/*
 * TrieNode.cpp
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

#include <iostream>
#include "TrieNode.h"
#include "Trie.h"

using namespace std;

TrieNode::TrieNode() {
    this->parent = nullptr;
    this->zeroBranch = nullptr;
    this->oneBranch = nullptr;
    this->wildcardBranch = nullptr;
    this->nextLevelTrie = nullptr;
    this->ruleSet = nullptr;
}

TrieNode::~TrieNode() {
//	fprintf(stdout, "Destroying TrieNode.\n");

    if (this->zeroBranch != nullptr) {
        // fprintf(stdout, "\tDestroying zeroBranch.\n");
        delete this->zeroBranch;
        this->zeroBranch = nullptr;
    }

    if (this->oneBranch != nullptr) {
        // fprintf(stdout, "\tDestroying oneBranch.\n");
        delete this->oneBranch;
        this->oneBranch = nullptr;
    }

    if (this->wildcardBranch != nullptr) {
        // fprintf(stdout, "\tDestroying wildcardBranch.\n");
        delete this->wildcardBranch;
        this->wildcardBranch = nullptr;
    }

    if (this->nextLevelTrie != nullptr) {
        // fprintf(stdout, "\tDestroying nextLevelTrie.\n");
        delete this->nextLevelTrie;
        this->nextLevelTrie = nullptr;
    }

    if (this->ruleSet != nullptr) {
        // fprintf(stdout, "\tDestroying ruleSet.\n");
        delete this->ruleSet;
        this->ruleSet = nullptr;
    }
}
