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
// This file is the header file for the MMTF parser for the C language.
//
// The authors of this code are: Julien Ferte (http://www.julienferte.com/),
// Anthony Bradley, Thomas Holder.
//
//
// Other contributors: Yana Valasatava, Gazal Kalyan, Alexander Rose.
//
// *************************************************************************

#ifndef MMTF_PARSER_H
#define MMTF_PARSER_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief MMTF spec version which this library implements
 */
#define MMTF_SPEC_VERSION_MAJOR 1
#define MMTF_SPEC_VERSION_MINOR 0

#define WITHCOUNT(name) \
    name; \
    size_t name##Count

/**
 * @brief Group (residue) level data store
 *
 * https://github.com/rcsb/mmtf/blob/HEAD/spec.md#group-data
 */
typedef struct {
    int32_t* formalChargeList;
    char** WITHCOUNT(atomNameList);
    char** WITHCOUNT(elementList);
    int32_t* WITHCOUNT(bondAtomList);
    int8_t* WITHCOUNT(bondOrderList);
    char* groupName;
    char singleLetterCode;
    char* chemCompType;
} MMTF_GroupType;

/**
 * @brief Entity type.
 *
 * https://github.com/rcsb/mmtf/blob/HEAD/spec.md#entitylist
 */
typedef struct {
    int32_t* WITHCOUNT(chainIndexList);
    char* description;
    char* type;
    char* sequence;
} MMTF_Entity;

/**
 * @brief Transformation definition for a set of chains.
 *
 * https://github.com/rcsb/mmtf/blob/HEAD/spec.md#bioassemblylist
 */
typedef struct {
    int32_t* WITHCOUNT(chainIndexList);
    float matrix[16];
} MMTF_Transform;

/**
 * @brief Data store for the biological assembly annotation.
 *
 * https://github.com/rcsb/mmtf/blob/HEAD/spec.md#bioassemblylist
 */
typedef struct {
    MMTF_Transform* WITHCOUNT(transformList);
    char* name;
} MMTF_BioAssembly;

/**
 * @brief Top level MMTF data container.
 *
 * https://github.com/rcsb/mmtf/blob/HEAD/spec.md#fields
 */
typedef struct {
    char* mmtfVersion;
    char* mmtfProducer;
    float unitCell[6];
    char* spaceGroup;
    char* structureId;
    char* title;
    char* depositionDate;
    char* releaseDate;
    MMTF_BioAssembly* WITHCOUNT(bioAssemblyList);
    MMTF_Entity* WITHCOUNT(entityList);
    char** WITHCOUNT(experimentalMethods);
    float resolution;
    float rFree;
    float rWork;
    int32_t numBonds;
    int32_t numAtoms;
    int32_t numGroups;
    int32_t numChains;
    int32_t numModels;
    MMTF_GroupType* WITHCOUNT(groupList);
    int32_t* WITHCOUNT(bondAtomList);
    int8_t* WITHCOUNT(bondOrderList);
    float* xCoordList;
    float* yCoordList;
    float* zCoordList;
    float* bFactorList;
    int32_t* atomIdList;
    char* altLocList;
    float* occupancyList;
    int32_t* groupIdList;
    int32_t* groupTypeList;
    int8_t* secStructList;
    char* insCodeList;
    int32_t* sequenceIndexList;
    char** WITHCOUNT(chainIdList);
    char** WITHCOUNT(chainNameList);
    int32_t* groupsPerChain;
    int32_t* chainsPerModel;
} MMTF_container;

/*
 * Macros for declaration generic initialization and destroying functions
 */

#define CODEGEN_MMTF_parser_TYPE(type) \
    type* type##_new(void); \
    void type##_init(type*); \
    void type##_destroy(type*); \
    void type##_clear(type*); \
    void type##_free(type*);

/**
 * @fn MMTF_container * MMTF_container_new()
 * @brief Allocate and initialize a new empty MMTF container
 * @return pointer which needs to be freed with ::MMTF_container_free
 *
 * @fn void MMTF_container_free(MMTF_container* container)
 * @brief Destroy the instance and free the memory
 * @param container A MMTF container which has been allocated with ::MMTF_container_new
 *
 * @fn void MMTF_container_init(MMTF_container * container)
 * @brief Initialize an empty MMTF container
 * @param container Uninitialized memory, needs to be destroyed with ::MMTF_container_destroy
 *
 * @fn void MMTF_container_destroy(MMTF_container* container)
 * @brief Destroy the instance and leave the memory uninitialized
 * @param container A MMTF container which has been initialized with ::MMTF_container_init
 */

// clang-format off
CODEGEN_MMTF_parser_TYPE(MMTF_container)
CODEGEN_MMTF_parser_TYPE(MMTF_BioAssembly)
CODEGEN_MMTF_parser_TYPE(MMTF_Transform)
CODEGEN_MMTF_parser_TYPE(MMTF_Entity)
CODEGEN_MMTF_parser_TYPE(MMTF_GroupType)

#undef CODEGEN_MMTF_parser_TYPE

/**
 * @brief Decode a MMTF_container from a string
 * @param[in] buffer file contents
 * @param[in] size file size
 * @param[out] container initialized but empty MMTF container to populate with data
 * @return true on success, false if an error occured
 */
bool MMTF_unpack_from_string(const char* buffer, size_t size, MMTF_container* container);
// clang-format on

/**
 * @brief Decode a MMTF_container from a file
 * @param[in] filename file path of an uncompressed .mmtf file
 * @param[out] container initialized but empty MMTF container to populate with data
 * @return true on success, false if an error occured
 */
bool MMTF_unpack_from_file(const char* filename, MMTF_container* container);

#undef WITHCOUNT

#ifdef __cplusplus
}
#endif
#endif

// vi:sw=4:expandtab
