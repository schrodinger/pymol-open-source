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
// This is a private header file for mmtf_parser.c and must not be included
// by other files.
//
// The authors of this code are: Julien Ferte (http://www.julienferte.com/),
// Anthony Bradley, Thomas Holder.
//
//
// Other contributors: Yana Valasatava, Alexander Rose.
//
// *************************************************************************

#ifndef MMTF_PRIVATE_H
#define MMTF_PRIVATE_H

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

#ifdef __cplusplus
extern "C" {
#endif

//*** Array converters
float* MMTF_parser_float_from_bytes(const char*, uint32_t, uint32_t*);
int8_t* MMTF_parser_int8_from_bytes(const char*, uint32_t, uint32_t*);
int16_t* MMTF_parser_int16_from_bytes(const char*, uint32_t, uint32_t*);
int32_t* MMTF_parser_int32_from_bytes(const char*, const uint32_t, uint32_t*);
char** MMTF_parser_strings_from_bytes(const char*, uint32_t, uint32_t, uint32_t*);

//*** Array decoders
int32_t* MMTF_parser_run_length_decode(const int32_t*, uint32_t, uint32_t*);
int32_t* MMTF_parser_delta_decode(const int32_t*, uint32_t, uint32_t*);
int32_t* MMTF_parser_recursive_indexing_decode_from_16(const int16_t*, uint32_t, uint32_t*);
int32_t* MMTF_parser_recursive_indexing_decode_from_8(const int8_t*, uint32_t, uint32_t*);
float* MMTF_parser_integer_decode_from_16(const int16_t*, uint32_t, int32_t, uint32_t*);
float* MMTF_parser_integer_decode_from_32(const int32_t*, uint32_t, int32_t, uint32_t*);

//*** Unpacking from MsgPack and applying strategy
char* MMTF_parser_fetch_string(const msgpack_object*);
int64_t MMTF_parser_fetch_int(const msgpack_object*);
float MMTF_parser_fetch_float(const msgpack_object*);

MMTF_Entity* MMTF_parser_fetch_entityList(const msgpack_object*, size_t*);

MMTF_GroupType* MMTF_parser_fetch_groupList(const msgpack_object*, size_t*);

MMTF_BioAssembly* MMTF_parser_fetch_bioAssemblyList(const msgpack_object*, size_t*);
MMTF_Transform* MMTF_parser_fetch_transformList(const msgpack_object*, size_t*);

#ifdef __cplusplus
}
#endif
#endif

// vi:sw=4:expandtab
