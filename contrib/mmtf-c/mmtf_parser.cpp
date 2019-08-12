// *************************************************************************
// Copyright [2016] [RCSB]
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//
// This file is mmtf_parser.c, holding the source code of the MMTF parser
//
// The authors of this code are: Julien Ferte (http://www.julienferte.com/),
// Anthony Bradley, Thomas Holder.
//
//
// Other contributors: Yana Valasatava, Gazal Kalyan, Alexander Rose.
//
// *************************************************************************

#define WIN32_LEAN_AND_MEAN
#define __STDC_LIMIT_MACROS

#include "mmtf_parser.h"

//*** Standard libs
#include <stdbool.h>

//*** MsgPack lib
#ifdef MMTF_MSGPACK_USE_CPP11
#include <msgpack.hpp>
#define msgpack_object msgpack::object
#define msgpack_object_kv msgpack::object_kv
#define msgpack_object_str msgpack::object_str
#define MMTF_MSGPACK_TYPE(T) msgpack::type::T
#else
#include <msgpack.h>
#define MMTF_MSGPACK_TYPE(T) MSGPACK_OBJECT_##T
#endif

#if MSGPACK_VERSION_MAJOR < 1
#error "msgpack-c >= 1.0 required (https://github.com/msgpack/msgpack-c)"
#else

//*** For the constant NAN
#include <math.h>

//*** Standard libs
#include <stdio.h>

// byteorder functions ("ntohl" etc.)
#ifdef WIN32
#include <Winsock2.h>
#else
#include <arpa/inet.h>
#endif

// Accept msgpack bin type for strings (emits warnings).
// Both types seem to have the same memory layout.
#define MMTF_ACCEPT_MSGPACK_BIN_AS_STR

// clang-format off
// typed array memory allocation
#define MALLOC_ARRAY(type, size) (type*) malloc((size) * sizeof(type))

/*
 * Type aliases for code generation
 */

#define TYPEALIAS_char      char
#define TYPEALIAS_int8      int8_t
#define TYPEALIAS_int32     int32_t
#define TYPEALIAS_float     float
#define TYPEALIAS_string    char*
#define TYPEALIAS_int       int
// clang-format on

enum {
    MMTF_TYPE_char,
    MMTF_TYPE_int8 = MMTF_TYPE_char,
    MMTF_TYPE_int16,
    MMTF_TYPE_int32,
    MMTF_TYPE_float,
    MMTF_TYPE_string
};

/*
 * Macros for null-pointer checking
 */

#define IF_NULL_PTRERROR_RETURN(ptr, returnvalue) \
    if (!ptr) { \
        fprintf(stderr, "Error in %s: NULL pointer.\n", __FUNCTION__); \
        return returnvalue; \
    }

#define IF_NULL_ALLOCERROR_RETURN(ptr, returnvalue) \
    if (!ptr) { \
        fprintf(stderr, "Error in %s: couldn't allocate memory.\n", __FUNCTION__); \
        return returnvalue; \
    }

#define IF_NOT_MULTIPLE_ERROR_RETURN(length, size, returnvalue) \
    if ((length) % (size) != 0) { \
        fprintf(stderr, "Error in %s: length %u is not a multiple of %u.\n", __FUNCTION__, length, size); \
        return returnvalue; \
    }

#define IF_NULL_ALLOCERROR_RETURN_NULL(ptr) \
    IF_NULL_ALLOCERROR_RETURN(ptr, NULL)

/*
 * Macros for iterating over a msgpack map
 */

#define MAP_ITERATE_BEGIN_RV(object, returnvalue) \
    if (object->type != MMTF_MSGPACK_TYPE(MAP)) { \
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not a map.\n", __FUNCTION__); \
        return returnvalue; \
    } \
    { \
        msgpack_object_kv* current_key_value = object->via.map.ptr; \
        msgpack_object_kv* last_key_value = current_key_value + object->via.map.size; \
        for (; current_key_value != last_key_value; ++current_key_value) { \
            const msgpack_object* key = &(current_key_value->key); \
            const msgpack_object* value = &(current_key_value->val); \
            if (key->type == MMTF_MSGPACK_TYPE(BIN)) { \
                fprintf(stderr, "Warning: map key of type BIN ('%.*s')\n", key->via.bin.size, key->via.bin.ptr); \
            } else if (key->type != MMTF_MSGPACK_TYPE(STR)) { \
                fprintf(stderr, "Warning: map key not of type str (type %d).\n", key->type); \
                continue; \
            }

#define MAP_ITERATE_BEGIN(object) \
    MAP_ITERATE_BEGIN_RV(object, )

#define MAP_ITERATE_END() \
    } \
    }

/*
 * Macros for inside the map iteration
 */

#define FETCH_AND_ASSIGN(this_, type, name) \
    if (MMTF_parser_compare_msgpack_string_char_array(&(key->via.str), #name)) { \
        this_->name = MMTF_parser_fetch_##type(value); \
        continue; \
    }

#define FETCH_AND_ASSIGN_DUMMYCOUNT(this_, type, name) \
    if (MMTF_parser_compare_msgpack_string_char_array(&(key->via.str), #name)) { \
        size_t _length; \
        this_->name = MMTF_parser_fetch_##type(value, &_length); \
        continue; \
    }

#define FETCH_AND_ASSIGN_WITHCOUNT(this_, type, name) \
    if (MMTF_parser_compare_msgpack_string_char_array(&(key->via.str), #name)) { \
        this_->name = MMTF_parser_fetch_##type(value, &(this_->name##Count)); \
        continue; \
    }

#define FETCH_AND_ASSIGN_ARRAY(this_, type, name) \
    if (MMTF_parser_compare_msgpack_string_char_array(&(key->via.str), #name)) { \
        size_t _length; \
        type* array = MMTF_parser_fetch_##type##_array(value, &_length); \
        if (array != NULL) { \
            size_t i; \
            for (i = 0; i < _length; ++i) { \
                this_->name[i] = array[i]; \
            } \
            free(array); \
        } \
    }

/*
 * Macros for generating generic initialization and destroying functions
 */

#define CODEGEN_MMTF_parser_TYPE_init(type) \
    void type##_init(type* result) { \
        memset(result, 0, sizeof(type)); \
    }

#define CODEGEN_MMTF_parser_TYPE_new(type) \
    type* type##_new(void) { \
        type* result = (type*)malloc(sizeof(type)); \
        IF_NULL_ALLOCERROR_RETURN_NULL(result); \
        type##_init(result); \
        return result; \
    }

#define CODEGEN_MMTF_parser_TYPE_clear(type) \
    void type##_clear(type* result) { \
        IF_NULL_PTRERROR_RETURN(result, ); \
        type##_destroy(result); \
        type##_init(result); \
    }

#define CODEGEN_MMTF_parser_TYPE_free(type) \
    void type##_free(type* thing) { \
        IF_NULL_PTRERROR_RETURN(thing, ); \
        type##_destroy(thing); \
        free(thing); \
    }

// clang-format off
#define CODEGEN_MMTF_parser_TYPE(type) \
    CODEGEN_MMTF_parser_TYPE_init(type) \
    CODEGEN_MMTF_parser_TYPE_new(type) \
    CODEGEN_MMTF_parser_TYPE_clear(type) \
    CODEGEN_MMTF_parser_TYPE_free(type)
// clang-format on

#define generic_destroy(ptr) \
    free(*(ptr))

#define FREE_LIST(type_, name) \
    if (name != NULL) { \
        size_t i; \
        for (i = 0; i < name##Count; ++i) { \
            type_##_destroy(name + i); \
        } \
        free(name); \
    }

/*
 * Macros for fetching and converting msgpack array objects
 */

#define CODEGEN_BODY_fetch_OBJECT_ARRAY(type_, RESULT_I_ASSIGN) \
    { \
        if (object->type != MMTF_MSGPACK_TYPE(ARRAY)) { \
            fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not an array.\n", __FUNCTION__); \
            return NULL; \
        } \
        const msgpack_object* iter = object->via.array.ptr; \
        (*length) = object->via.array.size; \
        const msgpack_object* iter_end = iter + (*length); \
        type_* result = MALLOC_ARRAY(type_, *length); \
        IF_NULL_ALLOCERROR_RETURN_NULL(result); \
        int i = 0; \
        for (; iter != iter_end; ++iter, ++i) { \
            RESULT_I_ASSIGN; \
        } \
        return result; \
    }

// clang-format off
#define CODEGEN_MMTF_parser_fetch_array(type_, RESULT_I_ASSIGN) \
    static TYPEALIAS_##type_* MMTF_parser_fetch_##type_##_array( \
            const msgpack_object* object, size_t* length) { \
        if (object->type == MMTF_MSGPACK_TYPE(BIN)) { \
            return (TYPEALIAS_##type_*) \
            MMTF_parser_fetch_typed_array(object, length, MMTF_TYPE_##type_); \
        } \
        CODEGEN_BODY_fetch_OBJECT_ARRAY(TYPEALIAS_##type_, RESULT_I_ASSIGN); \
    }

#define CODEGEN_MMTF_parser_fetch_List(type_, suffix) \
    static type_* MMTF_parser_fetch_##suffix##List( \
            const msgpack_object* object, size_t* length) { \
        CODEGEN_BODY_fetch_OBJECT_ARRAY(type_, { \
            MMTF_parser_put_##suffix(iter, result + i); \
        }) \
    }
// clang-format on

/*
 * Generate "initialize", "new", "empty" and "destroy" functions for MMTF struct types.
 */
// clang-format off
CODEGEN_MMTF_parser_TYPE(MMTF_container)
CODEGEN_MMTF_parser_TYPE(MMTF_BioAssembly)
CODEGEN_MMTF_parser_TYPE(MMTF_Transform)
CODEGEN_MMTF_parser_TYPE(MMTF_Entity)
CODEGEN_MMTF_parser_TYPE(MMTF_GroupType)

//*** Destroy the innner of a struct
void MMTF_container_destroy(MMTF_container* thing) {
    // clang-format on
    IF_NULL_PTRERROR_RETURN(thing, );

    FREE_LIST(MMTF_BioAssembly, thing->bioAssemblyList);
    FREE_LIST(MMTF_Entity, thing->entityList);
    FREE_LIST(generic, thing->experimentalMethods);
    FREE_LIST(MMTF_GroupType, thing->groupList);
    FREE_LIST(generic, thing->chainIdList);
    FREE_LIST(generic, thing->chainNameList);

    free(thing->mmtfVersion);
    free(thing->mmtfProducer);
    free(thing->spaceGroup);
    free(thing->structureId);
    free(thing->title);
    free(thing->depositionDate);
    free(thing->releaseDate);
    free(thing->bondAtomList);
    free(thing->bondOrderList);
    free(thing->xCoordList);
    free(thing->yCoordList);
    free(thing->zCoordList);
    free(thing->bFactorList);
    free(thing->atomIdList);
    free(thing->altLocList);
    free(thing->occupancyList);
    free(thing->groupIdList);
    free(thing->groupTypeList);
    free(thing->secStructList);
    free(thing->insCodeList);
    free(thing->sequenceIndexList);
    free(thing->groupsPerChain);
    free(thing->chainsPerModel);

    // PyMOL
    free(thing->pymolRepsList);
    free(thing->pymolColorList);
}
void MMTF_BioAssembly_destroy(MMTF_BioAssembly* bio_assembly) {
    IF_NULL_PTRERROR_RETURN(bio_assembly, );
    FREE_LIST(MMTF_Transform, bio_assembly->transformList);
    free(bio_assembly->name);
}
void MMTF_Transform_destroy(MMTF_Transform* transform) {
    IF_NULL_PTRERROR_RETURN(transform, );
    free(transform->chainIndexList);
}
void MMTF_Entity_destroy(MMTF_Entity* entity) {
    IF_NULL_PTRERROR_RETURN(entity, );
    free(entity->chainIndexList);
    free(entity->description);
    free(entity->type);
    free(entity->sequence);
}
void MMTF_GroupType_destroy(MMTF_GroupType* group_type) {
    IF_NULL_PTRERROR_RETURN(group_type, );
    FREE_LIST(generic, group_type->atomNameList);
    FREE_LIST(generic, group_type->elementList);
    free(group_type->formalChargeList);
    free(group_type->bondAtomList);
    free(group_type->bondOrderList);
    free(group_type->groupName);
    free(group_type->chemCompType);
}

//*** Array converters
// From bytes[] to float32[], int8[], int16[], int32[] and string

static inline
void assign_bigendian_4(void* dst, const char* src) {
    *((uint32_t*)dst) = ntohl(*((uint32_t*)src));
}

static inline
void assign_bigendian_2(void* dst, const char* src) {
    *((uint16_t*)dst) = ntohs(*((uint16_t*)src));
}

static
void array_copy_bigendian_4(void* dst, const char* src, size_t n) {
    size_t i;
    for (i = 0; i < n; i += 4) {
        assign_bigendian_4(((char*)dst) + i, src + i);
    }
}

static
void array_copy_bigendian_2(void* dst, const char* src, size_t n) {
    size_t i;
    for (i = 0; i < n; i += 2) {
        assign_bigendian_2(((char*)dst) + i, src + i);
    }
}

static
float* MMTF_parser_float_from_bytes(const char* input, uint32_t input_length, uint32_t* output_length) {
    IF_NOT_MULTIPLE_ERROR_RETURN(input_length, 4, NULL);

    (*output_length) = input_length / 4;

    float* output = MALLOC_ARRAY(float, *output_length);
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    array_copy_bigendian_4(output, input, input_length);

    return output;
}

static
int8_t* MMTF_parser_int8_from_bytes(const char* input, uint32_t input_length, uint32_t* output_length) {
    (*output_length) = input_length;

    int8_t* output = MALLOC_ARRAY(int8_t, *output_length);
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    memcpy(output, input, input_length);

    return output;
}

static
int16_t* MMTF_parser_int16_from_bytes(const char* input, uint32_t input_length, uint32_t* output_length) {
    IF_NOT_MULTIPLE_ERROR_RETURN(input_length, 2, NULL);

    (*output_length) = input_length / 2;

    int16_t* output = MALLOC_ARRAY(int16_t, (*output_length));
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    array_copy_bigendian_2(output, input, input_length);

    return output;
}

static
int32_t* MMTF_parser_int32_from_bytes(const char* input, const uint32_t input_length, uint32_t* output_length) {
    IF_NOT_MULTIPLE_ERROR_RETURN(input_length, 4, NULL);

    (*output_length) = input_length / 4;

    int32_t* output = MALLOC_ARRAY(int32_t, (*output_length));
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    array_copy_bigendian_4(output, input, input_length);

    return output;
}

static
char** MMTF_parser_strings_from_bytes(const char* input, uint32_t input_length, uint32_t parameter, uint32_t* output_length) {
    IF_NOT_MULTIPLE_ERROR_RETURN(input_length, parameter, NULL);

    (*output_length) = input_length / parameter;

    char** output = MALLOC_ARRAY(char*, (*output_length));
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    uint32_t i;
    for (i = 0; i < *output_length; ++i) {
        output[i] = MALLOC_ARRAY(char, parameter + 1);
        IF_NULL_ALLOCERROR_RETURN_NULL(output[i]);
        memcpy(output[i], input + (i * parameter), parameter);
        output[i][parameter] = 0;
    }

    return output;
}

//*** Array decoders
// Run-length decode
static
int32_t* MMTF_parser_run_length_decode(const int32_t* input, uint32_t input_length, uint32_t* output_length) {
    (*output_length) = 0;

    IF_NOT_MULTIPLE_ERROR_RETURN(input_length, 2, NULL);

    uint32_t i;
    int32_t value, number;
    for (i = 0; i < input_length; i += 2) {
        number = input[i + 1];

        (*output_length) += number;
    }

    int32_t* output = MALLOC_ARRAY(int32_t, (*output_length)); // The output needs to be freed by the calling process
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    int j = 0;
    int k;
    for (i = 0; i < input_length; i += 2) {
        value = input[i];
        number = input[i + 1];

        for (k = 0; k < number; ++k) {
            output[j] = value;
            ++j;
        }
    }

    return output;
}

// Delta decode
static
int32_t* MMTF_parser_delta_decode(const int32_t* input, uint32_t input_length, uint32_t* output_length) {
    (*output_length) = input_length;
    int32_t* output = MALLOC_ARRAY(int32_t, (*output_length)); // The output needs to be freed by the calling process
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    if (input_length > 0) {
        output[0] = input[0];
    }

    uint32_t i;
    for (i = 1; i < input_length; ++i) {
        output[i] = output[i - 1] + input[i];
    }

    return output;
}

// Recursive indexing decode
static
int32_t* MMTF_parser_recursive_indexing_decode_from_16(const int16_t* input, uint32_t input_length, uint32_t* output_length) {
    (*output_length) = 0;
    uint32_t i;
    for (i = 0; i < input_length; ++i) {

        if (input[i] != INT16_MAX && input[i] != INT16_MIN) {
            ++(*output_length);
        }
    }

    int32_t* output = (int32_t*)MALLOC_ARRAY(int32_t, (*output_length)); // The output needs to be freed by the calling process
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    size_t j = 0;

    if (input_length > 0) {
        output[0] = 0;
    }

    for (i = 0; i < input_length; ++i) {
        output[j] += input[i];

        if (input[i] != INT16_MAX && input[i] != INT16_MIN && j + 1 < *output_length) {
            ++j;
            output[j] = 0;
        }
    }

    return output;
}

static
int32_t* MMTF_parser_recursive_indexing_decode_from_8(const int8_t* input, uint32_t input_length, uint32_t* output_length) {
    (*output_length) = 0;
    uint32_t i;
    for (i = 0; i < input_length; ++i) {
        if (input[i] != INT8_MAX && input[i] != INT8_MIN) {
            ++(*output_length);
        }
    }

    int32_t* output = MALLOC_ARRAY(int32_t, (*output_length)); // The output needs to be freed by the calling process
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    size_t j = 0;
    output[j] = 0;

    for (i = 0; i < input_length; ++i) {
        output[j] += input[i];

        if (input[i] != INT8_MAX && input[i] != INT8_MIN && j + 1 < *output_length) {
            ++j;
            output[j] = 0;
        }
    }

    return output;
}

// Integer decoding
static
float* MMTF_parser_integer_decode_from_16(const int16_t* input, uint32_t input_length, int32_t parameter, uint32_t* output_length) {
    (*output_length) = input_length;
    float* output = (float*)MALLOC_ARRAY(float, (*output_length));
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    float parameter_float = (float)parameter;
    uint32_t i;
    for (i = 0; i < input_length; ++i) {
        output[i] = ((float)input[i]) / parameter_float;
    }

    return output;
}

static
float* MMTF_parser_integer_decode_from_32(const int32_t* input, uint32_t input_length, int32_t parameter, uint32_t* output_length) {
    (*output_length) = input_length;
    float* output = (float*)MALLOC_ARRAY(float, (*output_length));
    IF_NULL_ALLOCERROR_RETURN_NULL(output);

    float parameter_float = (float)parameter;
    uint32_t i;
    for (i = 0; i < input_length; ++i) {
        output[i] = ((float)input[i]) / parameter_float;
    }

    return output;
}

//*** Applying a decoding strategy for getting an array
static
void* MMTF_parser_decode_apply_strategy(const char* input,
        uint32_t input_length, uint32_t* output_length, int strategy,
        int32_t parameter, int* typecode) {
    switch (strategy) {
    case 1:
        *typecode = MMTF_TYPE_float;
        return MMTF_parser_float_from_bytes(input, input_length, output_length);
    case 2:
        *typecode = MMTF_TYPE_int8;
        return MMTF_parser_int8_from_bytes(input, input_length, output_length);
    case 3:
        *typecode = MMTF_TYPE_int16;
        return MMTF_parser_int16_from_bytes(input, input_length, output_length);
    case 4:
        *typecode = MMTF_TYPE_int32;
        return MMTF_parser_int32_from_bytes(input, input_length, output_length);
    case 5:
        *typecode = MMTF_TYPE_string;
        return MMTF_parser_strings_from_bytes(input, input_length, parameter, output_length);
    case 6: {
        // pass
    }
    case 7: {
        uint32_t step1_length;
        int32_t* step1 = MMTF_parser_int32_from_bytes(input, input_length, &step1_length);

        int32_t* output = MMTF_parser_run_length_decode(step1, step1_length, output_length);
        free(step1);

        if (strategy == 6) {
            uint32_t i = 0;
            char* char_output = MALLOC_ARRAY(char, (*output_length));
            IF_NULL_ALLOCERROR_RETURN_NULL(char_output);
            for (; i < *output_length; ++i) {
                char_output[i] = output[i];
            }
            free(output);

            *typecode = MMTF_TYPE_int8;
            return char_output;
        }

        *typecode = MMTF_TYPE_int32;
        return output;
    }
    case 8: {
        uint32_t step1_length;
        int32_t* step1 = MMTF_parser_int32_from_bytes(input, input_length, &step1_length);

        uint32_t step2_length;
        int32_t* step2 = MMTF_parser_run_length_decode(step1, step1_length, &step2_length);
        free(step1);

        int32_t* output = MMTF_parser_delta_decode(step2, step2_length, output_length);
        free(step2);

        *typecode = MMTF_TYPE_int32;
        return output;
    }
    case 9: {
        uint32_t step1_length;
        int32_t* step1 = MMTF_parser_int32_from_bytes(input, input_length, &step1_length);

        uint32_t step2_length;
        int32_t* step2 = MMTF_parser_run_length_decode(step1, step1_length, &step2_length);
        free(step1);

        float* output = MMTF_parser_integer_decode_from_32(step2, step2_length, parameter, output_length);
        free(step2);

        *typecode = MMTF_TYPE_float;
        return output;
    }
    case 10: {
        uint32_t step1_length;
        int16_t* step1 = MMTF_parser_int16_from_bytes(input, input_length, &step1_length);

        uint32_t step2_length;
        int32_t* step2 = MMTF_parser_recursive_indexing_decode_from_16(step1, step1_length, &step2_length);
        free(step1);

        uint32_t step3_length;
        int32_t* step3 = MMTF_parser_delta_decode(step2, step2_length, &step3_length);
        free(step2);

        float* output = MMTF_parser_integer_decode_from_32(step3, step3_length, parameter, output_length);
        free(step3);

        *typecode = MMTF_TYPE_float;
        return output;
    }
    case 11: {
        uint32_t step1_length;
        int16_t* step1 = MMTF_parser_int16_from_bytes(input, input_length, &step1_length);

        float* output = MMTF_parser_integer_decode_from_16(step1, step1_length, parameter, output_length);
        free(step1);

        *typecode = MMTF_TYPE_float;
        return output;
    }
    case 12: {
        uint32_t step1_length;
        int16_t* step1 = MMTF_parser_int16_from_bytes(input, input_length, &step1_length);

        uint32_t step2_length;
        int32_t* step2 = MMTF_parser_recursive_indexing_decode_from_16(step1, step1_length, &step2_length);
        free(step1);

        float* output = MMTF_parser_integer_decode_from_32(step2, step2_length, parameter, output_length);
        free(step2);

        *typecode = MMTF_TYPE_float;
        return output;
    }
    case 13: {
        uint32_t step1_length;
        int8_t* step1 = MMTF_parser_int8_from_bytes(input, input_length, &step1_length);

        uint32_t step2_length;
        int32_t* step2 = MMTF_parser_recursive_indexing_decode_from_8(step1, step1_length, &step2_length);
        free(step1);

        float* output = MMTF_parser_integer_decode_from_32(step2, step2_length, parameter, output_length);
        free(step2);

        *typecode = MMTF_TYPE_float;
        return output;
    }
    case 14: {
        uint32_t step1_length;
        int16_t* step1 = MMTF_parser_int16_from_bytes(input, input_length, &step1_length);

        int32_t* output = MMTF_parser_recursive_indexing_decode_from_16(step1, step1_length, output_length);
        free(step1);

        *typecode = MMTF_TYPE_int32;
        return output;
    }
    case 15: {
        uint32_t step1_length;
        int8_t* step1 = MMTF_parser_int8_from_bytes(input, input_length, &step1_length);

        int32_t* output = MMTF_parser_recursive_indexing_decode_from_8(step1, step1_length, output_length);
        free(step1);

        *typecode = MMTF_TYPE_int32;
        return output;
    }
    default: {
        fprintf(stderr, "Error in %s: %i does not refer to any strategy.\n", __FUNCTION__, strategy);
        return NULL;
    }
    }
}

/*
 * Copy string from 'object' to 'out'
 */
static
void MMTF_parser_put_string(const msgpack_object* object, char** out) {
    size_t string_size = object->via.str.size;
    char* result = (*out) = MALLOC_ARRAY(char, (string_size + 1));
    IF_NULL_ALLOCERROR_RETURN(result, );
    memcpy(result, object->via.str.ptr, string_size);
    result[string_size] = '\0';
}

//*** Unpacking from MsgPack and applying strategy
static
char* MMTF_parser_fetch_string(const msgpack_object* object) {
    switch (object->type) {
#ifdef MMTF_ACCEPT_MSGPACK_BIN_AS_STR
    case MMTF_MSGPACK_TYPE(BIN):
        fprintf(stderr, "Warning in %s: type BIN, expected STR ('%.*s')\n", __FUNCTION__,
                object->via.bin.size, object->via.bin.ptr);
#endif
    case MMTF_MSGPACK_TYPE(STR):
        break;
    default:
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not a string.\n", __FUNCTION__);
        return NULL;
    }

    char* result = NULL;
    MMTF_parser_put_string(object, &result);
    return result;
}

static
char MMTF_parser_fetch_char(const msgpack_object* object) {
    switch (object->type) {
#ifdef MMTF_ACCEPT_MSGPACK_BIN_AS_STR
    case MMTF_MSGPACK_TYPE(BIN):
        fprintf(stderr, "Warning in %s: type BIN, expected STR ('%.*s')\n", __FUNCTION__,
                object->via.bin.size, object->via.bin.ptr);
#endif
    case MMTF_MSGPACK_TYPE(STR):
        break;
    default:
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not a string.\n", __FUNCTION__);
        return '\0';
    }

    return *(object->via.str.ptr);
}

static
int64_t MMTF_parser_fetch_int(const msgpack_object* object) {
    int64_t result;

    if (object->type == MMTF_MSGPACK_TYPE(POSITIVE_INTEGER)) {
        result = object->via.u64;
    } else if (object->type == MMTF_MSGPACK_TYPE(NEGATIVE_INTEGER)) {
        result = object->via.i64;
    } else {
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not an integer.\n", __FUNCTION__);
        return 0;
    }

    return result;
}

static
float MMTF_parser_fetch_float(const msgpack_object* object) {
    switch (object->type) {
    case /* FLOAT64 */ MMTF_MSGPACK_TYPE(FLOAT):
#if MSGPACK_VERSION_MAJOR >= 2
    case /* FLOAT32 */ 0x0a: // msgpack-c >= 2.1
#endif
        return (float)object->via.f64;
    default:
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not a float.\n", __FUNCTION__);
        return NAN;
    }
}

/*
 * Fetch a compressed typed array
 */
static
void* MMTF_parser_fetch_typed_array(const msgpack_object* object, size_t* length, int typecode) {
    if (object->type != MMTF_MSGPACK_TYPE(BIN)) {
        fprintf(stderr, "Error in %s: the entry encoded in the MMTF is not binary data.\n", __FUNCTION__);
        return NULL;
    }

    const char* bytes = object->via.bin.ptr;

    int32_t strategy, len_int32, parameter;
    assign_bigendian_4(&strategy, bytes);
    assign_bigendian_4(&len_int32, bytes + 4);
    assign_bigendian_4(&parameter, bytes + 8);

    *length = len_int32;

    uint32_t out_length;
    int typecheck;
    void* result = MMTF_parser_decode_apply_strategy(bytes + 12,
            object->via.bin.size - 12, &out_length, strategy, parameter,
            &typecheck);

    if (typecode != typecheck) {
        fprintf(stderr, "Error in %s: typecode mismatch %d %d\n",
                __FUNCTION__, typecode, typecheck);
        return NULL;
    }

    if (out_length != *length) {
        fprintf(stderr, "Error in %s: length mismatch %u %u\n",
                __FUNCTION__, out_length, (unsigned)*length);
        return NULL;
    }

    return result;
}

/*
 * Fetch a typed array.
 */
// clang-format off
CODEGEN_MMTF_parser_fetch_array(char,   result[i] = iter->via.u64)
CODEGEN_MMTF_parser_fetch_array(int8,   result[i] = iter->via.u64)
CODEGEN_MMTF_parser_fetch_array(int32,  result[i] = iter->via.u64)
CODEGEN_MMTF_parser_fetch_array(float,  result[i] = iter->via.f64)
CODEGEN_MMTF_parser_fetch_array(string, MMTF_parser_put_string(iter, result + i))

static
bool MMTF_parser_compare_msgpack_string_char_array(const msgpack_object_str* m_string, const char* string) {
    // clang-format on
    return (m_string->size == strlen(string) && strncmp(m_string->ptr, string, m_string->size) == 0);
}

static
void MMTF_parser_put_entity(const msgpack_object* object, MMTF_Entity* entity) {
    MMTF_Entity_init(entity);
    MAP_ITERATE_BEGIN(object);
    FETCH_AND_ASSIGN(entity, string, description);
    FETCH_AND_ASSIGN(entity, string, type);
    FETCH_AND_ASSIGN_WITHCOUNT(entity, int32_array, chainIndexList);
    FETCH_AND_ASSIGN(entity, string, sequence);
    MAP_ITERATE_END();
}

static
void MMTF_parser_put_group(const msgpack_object* object, MMTF_GroupType* group_type) {
    MMTF_GroupType_init(group_type);
    MAP_ITERATE_BEGIN(object);
    FETCH_AND_ASSIGN_DUMMYCOUNT(group_type, int32_array, formalChargeList);
    FETCH_AND_ASSIGN_WITHCOUNT(group_type, string_array, atomNameList);
    FETCH_AND_ASSIGN_WITHCOUNT(group_type, string_array, elementList);
    FETCH_AND_ASSIGN_WITHCOUNT(group_type, int32_array, bondAtomList);
    FETCH_AND_ASSIGN_WITHCOUNT(group_type, int8_array, bondOrderList);
    FETCH_AND_ASSIGN(group_type, string, groupName);
    FETCH_AND_ASSIGN(group_type, char, singleLetterCode);
    FETCH_AND_ASSIGN(group_type, string, chemCompType);
    MAP_ITERATE_END();
}

static
MMTF_Transform* MMTF_parser_fetch_transformList(const msgpack_object*, size_t*);

static
void MMTF_parser_put_bioAssembly(const msgpack_object* object, MMTF_BioAssembly* bio_assembly) {
    MMTF_BioAssembly_init(bio_assembly);
    MAP_ITERATE_BEGIN(object);
    FETCH_AND_ASSIGN(bio_assembly, string, name);
    FETCH_AND_ASSIGN_WITHCOUNT(bio_assembly, transformList, transformList);
    MAP_ITERATE_END();
}

static
void MMTF_parser_put_transform(const msgpack_object* object, MMTF_Transform* transform) {
    MAP_ITERATE_BEGIN(object);
    FETCH_AND_ASSIGN_WITHCOUNT(transform, int32_array, chainIndexList);
    FETCH_AND_ASSIGN_ARRAY(transform, float, matrix);
    MAP_ITERATE_END();
}

// clang-format off
CODEGEN_MMTF_parser_fetch_List(MMTF_Entity, entity)
CODEGEN_MMTF_parser_fetch_List(MMTF_GroupType, group)
CODEGEN_MMTF_parser_fetch_List(MMTF_BioAssembly, bioAssembly)
CODEGEN_MMTF_parser_fetch_List(MMTF_Transform, transform)

static
bool MMTF_unpack_from_msgpack_object(const msgpack_object* object, MMTF_container* thing) {
    // clang-format on
    int version_major;

    MAP_ITERATE_BEGIN_RV(object, false);
    FETCH_AND_ASSIGN(thing, string, mmtfVersion);
    MAP_ITERATE_END();

    // check for MMTF major version (semantic versioning)
    if (thing->mmtfVersion &&
            sscanf(thing->mmtfVersion, "%d", &version_major) == 1 &&
            version_major > MMTF_SPEC_VERSION_MAJOR) {
        fprintf(stderr, "Error: unsupported MMTF version '%s'.\n", thing->mmtfVersion);
        return false;
    }

    MAP_ITERATE_BEGIN_RV(object, false);
    FETCH_AND_ASSIGN(thing, string, mmtfProducer);
    FETCH_AND_ASSIGN(thing, string, spaceGroup);
    FETCH_AND_ASSIGN(thing, string, structureId);
    FETCH_AND_ASSIGN(thing, string, title);
    FETCH_AND_ASSIGN(thing, string, depositionDate);
    FETCH_AND_ASSIGN(thing, string, releaseDate);
    FETCH_AND_ASSIGN(thing, int, numBonds);
    FETCH_AND_ASSIGN(thing, int, numAtoms);
    FETCH_AND_ASSIGN(thing, int, numGroups);
    FETCH_AND_ASSIGN(thing, int, numChains);
    FETCH_AND_ASSIGN(thing, int, numModels);
    FETCH_AND_ASSIGN(thing, float, resolution);
    FETCH_AND_ASSIGN(thing, float, rFree);
    FETCH_AND_ASSIGN(thing, float, rWork);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, entityList, entityList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, bioAssemblyList, bioAssemblyList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, groupList, groupList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, int32_array, bondAtomList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, int8_array, bondOrderList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, string_array, chainIdList);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, string_array, chainNameList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, groupTypeList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, groupIdList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, sequenceIndexList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, atomIdList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, char_array, insCodeList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, char_array, altLocList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int8_array, secStructList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, float_array, bFactorList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, float_array, xCoordList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, float_array, yCoordList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, float_array, zCoordList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, float_array, occupancyList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, chainsPerModel);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, groupsPerChain);
    FETCH_AND_ASSIGN_WITHCOUNT(thing, string_array, experimentalMethods);
    FETCH_AND_ASSIGN_ARRAY(thing, float, unitCell);

    // PyMOL
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, pymolRepsList);
    FETCH_AND_ASSIGN_DUMMYCOUNT(thing, int32_array, pymolColorList);

    MAP_ITERATE_END();

    return true;
}

/*
 * Decode a MMTF_container from a string
 */
bool MMTF_unpack_from_string(const char* buffer, size_t msgsize, MMTF_container* thing) {
    bool status;

#ifdef MMTF_MSGPACK_USE_CPP11

    auto oh = msgpack::unpack(buffer, msgsize);

    status = MMTF_unpack_from_msgpack_object(&oh.get(), thing);

#else

    msgpack_zone mempool;
    msgpack_zone_init(&mempool, 2048);
    msgpack_object deserialized;
    msgpack_unpack(buffer, msgsize, NULL, &mempool, &deserialized);

    status = MMTF_unpack_from_msgpack_object(&deserialized, thing);

    msgpack_zone_destroy(&mempool);

#endif

    return status;
}

/*
 * Decode a MMTF container from a file
 */
bool MMTF_unpack_from_file(const char* name, MMTF_container* thing) {
    FILE* file;
    char* buffer;
    size_t fileLen;
    bool status;

    //*** Open file
    file = fopen(name, "rb");
    if (!file) {
        fprintf(stderr, "Error in %s: unable to open file %s.\n", __FUNCTION__, name);
        return false;
    }

    //*** Get file length
    fseek(file, 0, SEEK_END);
    fileLen = ftell(file);
    fseek(file, 0, SEEK_SET);

    //*** Allocate memory
    buffer = (char*)malloc(fileLen + 1);
    if (!buffer) {
        fprintf(stderr, "Error in %s: couldn't allocate memory.\n", __FUNCTION__);
        fclose(file);
        return false;
    }

    //*** Read file contents into buffer
    fread(buffer, fileLen, 1, file);
    fclose(file);

    status = MMTF_unpack_from_string(buffer, fileLen, thing);

    free(buffer);

    return status;
}

#endif // MSGPACK_VERSION_MAJOR check

// vi:sw=4:expandtab
