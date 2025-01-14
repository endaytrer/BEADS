/*
 * EquivalenceClass.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: khurshid
 *
 * This file is protected by the VeriFlow Research License Agreement
 * available at http://www.cs.illinois.edu/~khurshi1/projects/veriflow/veriflow-research-license-agreement.txt.
 * A copy of this agreement is also included in this package.
 *
 * Copyright (c) 2012-2013 by
 * The Board of Trustees of the University of Illinois.
 * All rights reserved.
 */

#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include "EquivalenceClass.h"
#include "VeriFlow.h"

using namespace std;

EquivalenceClass::EquivalenceClass() {
    this->clear();
}

EquivalenceClass::EquivalenceClass(const uint64_t *lb, const uint64_t *ub) {
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        if (lb[i] > ub[i]) {
            fprintf(stderr,
                    "[EquivalenceClass::EquivalenceClass] Invalid boundary values for field %d (lb > ub). Terminating process.\n",
                    i);
            fprintf(stderr, "\tLower bound: %lu\n", lb[i]);
            fprintf(stderr, "\tUpper bound: %lu\n", ub[i]);

            exit(1);
        }

        this->lowerBound[i] = lb[i];
        this->upperBound[i] = ub[i];
    }
}

uint64_t EquivalenceClass::getRange(int fieldIndex) const {
    return (this->upperBound[fieldIndex] - this->lowerBound[fieldIndex]);
}

bool EquivalenceClass::equals(const EquivalenceClass &other) const {
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        if ((this->lowerBound[i] != other.lowerBound[i]) || (this->upperBound[i] != other.upperBound[i])) {
            return false;
        }
    }

    return true;
}

bool EquivalenceClass::operator==(const EquivalenceClass &other) const {
    return this->equals(other);
}
//NOTE: other is the subset of *this
bool EquivalenceClass::subsumes(const EquivalenceClass &other) const {
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        if ((this->lowerBound[i] <= other.lowerBound[i] && this->upperBound[i] >= other.upperBound[i])) {
            //Okay
        } else {
            return false;
        }
    }

    return true;

}

int EquivalenceClass::operator()() const { //INFO: for hashing.
    int retVal = 0;
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        retVal += (this->lowerBound[i] + this->upperBound[i]);
    }

    return retVal;
}

void EquivalenceClass::clear() {
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        this->lowerBound[i] = 0;
        this->upperBound[i] = 0;
    }
}

string EquivalenceClass::toString() const {
    char buffer[1024];
    sprintf(buffer, "nw_src (%s-%s), nw_dst (%s-%s)",
            ::getIpValueAsString(this->lowerBound[NW_SRC]).c_str(),
            ::getIpValueAsString(this->upperBound[NW_SRC]).c_str(),
            ::getIpValueAsString(this->lowerBound[NW_DST]).c_str(),
            ::getIpValueAsString(this->upperBound[NW_DST]).c_str());

    string retVal = buffer;
    retVal += ", ";
    // nw-proto
    sprintf(buffer, "nw_proto(%lld-%lld)", this->lowerBound[NW_PROTO], this->upperBound[NW_PROTO]);
    retVal += buffer;
    retVal += ", ";


    // tp_src, tp_dst

    sprintf(buffer, "tp_src(%lld-%lld), tp_dst(%lld-%lld)", this->lowerBound[TP_SRC], this->upperBound[TP_SRC], this->lowerBound[TP_DST], this->upperBound[TP_DST]);
    retVal += buffer;
    return retVal;
}

uint64_t EquivalenceClass::getMaxValue(FieldIndex index) {
    switch (index) {
        case IN_PORT:
            return 0xFFFF;

        case DL_SRC:
            return 0xFFFFFFFFFFFFLLU;

        case DL_DST:
            return 0xFFFFFFFFFFFFLLU;

        case DL_TYPE:
            return 0xFFFF;

        case DL_VLAN:
            return 0xFFF;

        case DL_VLAN_PCP:
            return 0x7;

        case MPLS_LABEL:
            return 0xFFFFF;

        case MPLS_TC:
            return 0x7;

        case NW_SRC:
            return 0xFFFFFFFF;

        case NW_DST:
            return 0xFFFFFFFF;

        case NW_PROTO:
            return 0xFF;

        case NW_TOS:
            return 0x3F;

        case TP_SRC:
            return 0xFFFF;

        case TP_DST:
            return 0xFFFF;

        case METADATA:
            return 0xFFFFFFFFFFFFFFFFLLU;

        default:
            fprintf(stderr, "[EquivalenceClass::getMaxValue] Wrong field index (%d) specified. Terminating process.\n",
                    index);
            exit(1);
            return 0x0; // Something went wrong; wrong field index was specified. Unreachable code due to the presence of exit(1).
            break;
    }
}
