/*
 * Trie.cpp
 *
 *  Created on: Jul 4, 2012
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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stack>
#include "Trie.h"
#include "EquivalenceClass.h"
#include "VeriFlow.h"

using namespace std;

extern int mode;
extern vector<string> endhosts;

Trie::Trie(FieldIndex fi) {
    if (fi >= ALL_FIELD_INDEX_END_MARKER) {
        fprintf(stderr, "[Trie::Trie] Invalid field index (%d). Terminating process.\n", fi);
        exit(1);
    }

    this->root = nullptr;
    this->fieldIndex = fi;
    this->totalRuleCount = 0;
}

Trie::~Trie() {
    if (this->root != nullptr) {
        delete root;
        this->root = nullptr;
        this->totalRuleCount = 0;
    }
}

FieldIndex Trie::getFieldIndex() {
    return this->fieldIndex;
}

unsigned int Trie::getFieldWidth() {
    return ::fieldWidth[this->fieldIndex];
}

uint64_t Trie::getIntValue(FieldIndex index, const string &valueOrMask) {
    if ((index == DL_SRC) || (index == DL_DST)) {
        return ::getMacValueAsInt(valueOrMask);
    } else if ((index == NW_SRC) || (index == NW_DST)) {
        return ::getIpValueAsInt(valueOrMask);
    } else {
        return strtoul(valueOrMask.c_str(), nullptr, 10);

        // fprintf(stderr, "[Trie::getIntValue] Invalid field index (%d). Terminating process.\n", index);
        // exit(1);
    }
}

TrieNode *Trie::findNode(const string &fieldValue, const string &fieldMask) {
    if ((this->root == nullptr) || (this->totalRuleCount == 0)) {
        return nullptr;
    }

    uint64_t fieldValueInt = Trie::getIntValue(this->fieldIndex, fieldValue);
    uint64_t fieldMaskInt = Trie::getIntValue(this->fieldIndex, fieldMask);
    uint64_t maskedFieldValue = fieldValueInt & fieldMaskInt;

    TrieNode *currentNode = this->root;

    for (unsigned int i = 0; i < ::fieldWidth[this->fieldIndex]; i++) {
        uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[this->fieldIndex] - 1) - i);
        if ((fieldMaskInt & maskBit) == 0) // wildcard bit
        {
            if (currentNode->wildcardBranch == nullptr) {
                return nullptr;
            }

            currentNode = currentNode->wildcardBranch;
        } else {
            if ((maskedFieldValue & maskBit) == 0) // zero bit
            {
                if (currentNode->zeroBranch == nullptr) {
                    return nullptr;
                }

                currentNode = currentNode->zeroBranch;
            } else // one bit
            {
                if (currentNode->oneBranch == nullptr) {
                    return nullptr;
                }

                currentNode = currentNode->oneBranch;
            }
        }
    }

    return currentNode;
}

TrieNode *Trie::addRule(const Rule &rule) {
    if (this->root == nullptr) {
        this->root = new TrieNode;
        if (this->root == nullptr) {
            fprintf(stderr, "[Trie::addRule] Memory allocation error (this->root == nullptr). Terminating process.\n");
            exit(1);
        }
    }

    uint64_t fieldValue = Trie::getIntValue(this->fieldIndex, rule.fieldValue[this->fieldIndex]);
    uint64_t fieldMask = Trie::getIntValue(this->fieldIndex, rule.fieldMask[this->fieldIndex]);
    uint64_t maskedFieldValue = fieldValue & fieldMask;

    TrieNode *currentNode = this->root;

    for (unsigned int i = 0; i < ::fieldWidth[this->fieldIndex]; i++) {
        uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[this->fieldIndex] - 1) - i);
        if ((fieldMask & maskBit) == 0) // wildcard bit
        {
            if (currentNode->wildcardBranch == nullptr) {
                currentNode->wildcardBranch = new TrieNode;
                if (currentNode->wildcardBranch == nullptr) {
                    fprintf(stderr,
                            "[Trie::addRule] Memory allocation error (currentNode->wildcardBranch == nullptr). Terminating process.\n");
                    exit(1);
                }

                currentNode->wildcardBranch->parent = currentNode;
            }

            currentNode = currentNode->wildcardBranch;
        } else {
            if ((maskedFieldValue & maskBit) == 0) // zero bit
            {
                if (currentNode->zeroBranch == nullptr) {
                    currentNode->zeroBranch = new TrieNode;
                    if (currentNode->zeroBranch == nullptr) {
                        fprintf(stderr,
                                "[Trie::addRule] Memory allocation error (currentNode->zeroBranch == nullptr). Terminating process.\n");
                        exit(1);
                    }

                    currentNode->zeroBranch->parent = currentNode;
                }

                currentNode = currentNode->zeroBranch;
            } else // one bit
            {
                if (currentNode->oneBranch == nullptr) {
                    currentNode->oneBranch = new TrieNode;
                    if (currentNode->oneBranch == nullptr) {
                        fprintf(stderr,
                                "[Trie::addRule] Memory allocation error (currentNode->oneBranch == nullptr). Terminating process.\n");
                        exit(1);
                    }

                    currentNode->oneBranch->parent = currentNode;
                }

                currentNode = currentNode->oneBranch;
            }
        }
    }

    return currentNode;
}

void Trie::removeRule(TrieNode *node) { // since rules can only be added to leaf nodes, if parent have more than one child, delete parent as well.
    TrieNode *parent = node->parent;
    while (parent != nullptr) {
        if (
                ((parent->zeroBranch == node)
                 && (parent->oneBranch == nullptr)
                 && (parent->wildcardBranch == nullptr))
                ||
                ((parent->oneBranch == node)
                 && (parent->zeroBranch == nullptr)
                 && (parent->wildcardBranch == nullptr))
                ||
                ((parent->wildcardBranch == node)
                 && (parent->zeroBranch == nullptr)
                 && (parent->oneBranch == nullptr))
                ) {
            // fprintf(stdout, "Climbing up one level\n");
            TrieNode *tempParent = parent->parent;
            node = parent;
            parent = tempParent;
        } else {
            if (parent->zeroBranch == node) {
                // fprintf(stdout, "Deleting zero branch...\n");
                parent->zeroBranch = nullptr;
            } else if (parent->oneBranch == node) {
                // fprintf(stdout, "Deleting one branch...\n");
                parent->oneBranch = nullptr;
            } else if (parent->wildcardBranch == node) {
                // fprintf(stdout, "Deleting wildcard branch...\n");
                parent->wildcardBranch = nullptr;
            }

            delete node; // will cause chain reaction

            break;
        }
    }

    --(this->totalRuleCount);

    if (parent == nullptr) {
        delete this->root;
        // fprintf(stdout, "Destroyed root.\n");
        this->root = nullptr;
        this->totalRuleCount = 0;
    }
}

// Allow at most 1 occurrence of each value in each list.
void Trie::addToBoundList(uint64_t lowerBound, uint64_t upperBound) {
    if (this->lowerBoundMap.find(lowerBound) == this->lowerBoundMap.end()) {
        this->lowerBoundList.push_front(lowerBound);
        this->lowerBoundMap[lowerBound] = upperBound;
    } else {
        if (this->lowerBoundMap[lowerBound] < upperBound) {
            this->lowerBoundMap[lowerBound] = upperBound;
        }
    }

    if (this->upperBoundMap.find(upperBound) == this->upperBoundMap.end()) {
        this->upperBoundList.push_front(upperBound);
        this->upperBoundMap[upperBound] = lowerBound;
    } else {
        if (this->upperBoundMap[upperBound] > lowerBound) {
            this->upperBoundMap[upperBound] = lowerBound;
        }
    }
}

// Allow at most 1 occurrence of each value in each list.
void Trie::addToBoundList(uint64_t lowerBound, uint64_t upperBound,
                          list<uint64_t> &lowerBoundList, list<uint64_t> &upperBoundList,
                          unordered_map<uint64_t, uint64_t> &lowerBoundMap,
                          unordered_map<uint64_t, uint64_t> &upperBoundMap) {
    if (lowerBoundMap.find(lowerBound) == lowerBoundMap.end()) {
        lowerBoundList.push_front(lowerBound);
        lowerBoundMap[lowerBound] = upperBound;
    } else {
        if (lowerBoundMap[lowerBound] < upperBound) {
            lowerBoundMap[lowerBound] = upperBound;
        }
    }

    if (upperBoundMap.find(upperBound) == upperBoundMap.end()) {
        upperBoundList.push_front(upperBound);
        upperBoundMap[upperBound] = lowerBound;
    } else {
        if (upperBoundMap[upperBound] > lowerBound) {
            upperBoundMap[upperBound] = lowerBound;
        }
    }
}
EquivalenceRange Trie::getEquivalenceRangeAndTraverse(const Rule &rule, TrieNode *&out_Node) {
    // Traverse until wildcard appears.
    uint64_t fieldValue = Trie::getIntValue(this->fieldIndex, rule.fieldValue[this->fieldIndex]);
    uint64_t fieldMask = Trie::getIntValue(this->fieldIndex, rule.fieldMask[this->fieldIndex]);
    uint64_t maskedFieldValue = fieldValue & fieldMask;

    out_Node = this->root;
    EquivalenceRange range;
    for (unsigned int i = 0; i < ::fieldWidth[this->fieldIndex]; i++) {
        uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[this->fieldIndex] - 1) - i);
        if ((fieldMask & maskBit) != 0) // non-wildcard bit
        {
            if ((maskedFieldValue & maskBit) == 0) // zero bit
            {
                if (out_Node->zeroBranch != nullptr) {
                    out_Node = out_Node->zeroBranch;
                } else {
                    fprintf(stderr,
                            "[Trie::getEquivalenceClasses] Error: out_Node->zeroBranch cannot be nullptr. Terminating process.\n");
                    exit(1);
                }
            } else // one bit
            {
                range.lowerBound |= 1;
                range.upperBound |= 1;

                if (out_Node->oneBranch != nullptr) {
                    out_Node = out_Node->oneBranch;
                } else {
                    fprintf(stderr,
                            "[Trie::getEquivalenceClasses] Error: out_Node->oneBranch cannot be nullptr. Terminating process.\n");
                    exit(1);
                }
            }

            range.lowerBound <<= 1;
            range.upperBound <<= 1;
        } else // wildcard bit
        {
            break;
        }
    }
    return range;
}
void Trie::getEquivalenceClasses(const Rule &rule, vector<EquivalenceClass> &vPacketClasses) {
    vPacketClasses.clear();
    if (this->root == nullptr) {
        return;
    }

    TrieNode *currentNode;
    EquivalenceRange range = getEquivalenceRangeAndTraverse(rule, currentNode);
    this->lowerBoundList.clear();
    this->lowerBoundMap.clear();

    this->upperBoundList.clear();
    this->upperBoundMap.clear();

    // Perform DFS (Depth-first Search)

    stack<TrieNode *> nodes;
    nodes.push(currentNode);
    stack<EquivalenceRange> ranges;
    ranges.push(range);

    while (!nodes.empty()) {
        // Since we only reach the state of first wildcard, we have to count all nodes below.
        // DFS and add all leaves to ranges.
        currentNode = nodes.top();
        nodes.pop();
        EquivalenceRange tempRange = ranges.top();
        ranges.pop();

        // fprintf(stdout, "Exploring node...\n");

        if (currentNode == nullptr) {
            fprintf(stderr, "[Trie::getEquivalenceClasses] Invalid node (node = nullptr) found. Terminating process.\n");
            exit(1);
        }

        if (currentNode->oneBranch != nullptr) {
            EquivalenceRange oneRange = tempRange;
            oneRange.lowerBound |= 1;
            oneRange.upperBound |= 1;
            oneRange.lowerBound <<= 1;
            oneRange.upperBound <<= 1;
            ranges.push(oneRange);

            // fprintf(stdout, "Pushing oneBranch...\n");
            nodes.push(currentNode->oneBranch);
        }

        if (currentNode->zeroBranch != nullptr) {
            EquivalenceRange zeroRange = tempRange;
            // zeroRange.lowerBound |= 0;
            // zeroRange.upperBound |= 0;
            zeroRange.lowerBound <<= 1;
            zeroRange.upperBound <<= 1;
            ranges.push(zeroRange);

            // fprintf(stdout, "Pushing zeroBranch...\n");
            nodes.push(currentNode->zeroBranch);
        }

        if (currentNode->wildcardBranch != nullptr) {
            EquivalenceRange wildcardRange = tempRange;
            // wildcardRange.lowerBound |= 0;
            wildcardRange.upperBound |= 1;
            wildcardRange.lowerBound <<= 1;
            wildcardRange.upperBound <<= 1;
            ranges.push(wildcardRange);

            // fprintf(stdout, "Pushing wildcardBranch...\n");
            nodes.push(currentNode->wildcardBranch);
        }

        if ((currentNode->ruleSet != nullptr) || (currentNode->nextLevelTrie != nullptr)) {
            // Reached a leaf.

            tempRange.lowerBound >>= 1;
            tempRange.upperBound >>= 1;

            this->addToBoundList(tempRange.lowerBound, tempRange.upperBound);
        } else {
            // Do nothing.
        }
    }

    this->lowerBoundList.sort();
    this->upperBoundList.sort();

    // fprintf(stdout, "Generating equivalent packets...\n");

    uint64_t lowerBound = 0;
    uint64_t upperBound = 0;

    if (!this->lowerBoundList.empty()) {
        lowerBound = upperBound = this->lowerBoundList.front();
    }

    while (!this->upperBoundList.empty()) {
        if (!this->lowerBoundList.empty()) {
            lowerBound = this->lowerBoundList.front();
            if ((lowerBound > (upperBound + 1))
                /* && (this->upperBoundList.front() > (upperBound + 1))
                && (this->upperBoundMap[this->upperBoundList.front()] < (upperBound + 1)) */) {
                lowerBound = upperBound + 1;

                if (this->lowerBoundList.front() <= this->upperBoundList.front()) {
                    upperBound = this->lowerBoundList.front() - 1;
                } else {
                    upperBound = this->upperBoundList.front();
                    this->upperBoundList.pop_front();
                }

                EquivalenceClass packetClass;
                packetClass.lowerBound[this->fieldIndex] = lowerBound;
                packetClass.upperBound[this->fieldIndex] = upperBound;

                vPacketClasses.push_back(packetClass);

                continue;
            } else {
                this->lowerBoundList.pop_front();
            }
        } else {
            lowerBound = upperBound + 1;
        }

        if (!this->lowerBoundList.empty()) {
            if (this->lowerBoundList.front() <= this->upperBoundList.front()) {
                upperBound = this->lowerBoundList.front() - 1;
            } else {
                upperBound = this->upperBoundList.front();
                this->upperBoundList.pop_front();
            }
        } else {
            upperBound = this->upperBoundList.front();
            this->upperBoundList.pop_front();
        }

        EquivalenceClass packetClass;
        packetClass.lowerBound[this->fieldIndex] = lowerBound;
        packetClass.upperBound[this->fieldIndex] = upperBound;

        vPacketClasses.push_back(packetClass);
    }
}

void Trie::getNextLevelEquivalenceClasses(FieldIndex currentFieldIndex, uint64_t lb, const Rule &rule,
                                          const vector<Trie *> &vInputTries, vector<EquivalenceClass> &vPacketClasses,
                                          vector<Trie *> &vOutputTries) {
    vPacketClasses.clear();
    vOutputTries.clear();

    FieldIndex nextFieldIndex = (FieldIndex) (currentFieldIndex + 1);
    if (nextFieldIndex >= ALL_FIELD_INDEX_END_MARKER) {
        fprintf(stderr,
                "[Trie::getNextLevelEquivalenceClasses] Method called on wrong trie (field index = %d). Terminating process.\n",
                currentFieldIndex);
        exit(1);
    }

    uint64_t fieldValue = lb;
    uint64_t fieldMask = EquivalenceClass::getMaxValue(currentFieldIndex);

    uint64_t maskedFieldValue = fieldValue & fieldMask;

    for (auto inputTrie : vInputTries) {
        if (inputTrie->getTotalRuleCount() == 0) {
            continue;
        }

        vector<TrieNode *> vCurrentLevelNodes, vNextLevelNodes; // Layer order traversal
        vCurrentLevelNodes.push_back(inputTrie->root);
        TrieNode *currentNode;

        for (unsigned int i = 0; i < ::fieldWidth[currentFieldIndex]; i++) {
            while (!vCurrentLevelNodes.empty()) {
                currentNode = vCurrentLevelNodes.at(0);
                vCurrentLevelNodes.erase(vCurrentLevelNodes.begin());

                if (currentNode == nullptr) {
                    fprintf(stderr,
                            "[Trie::getNextLevelEquivalenceClasses] Invalid node (node = nullptr) found. Terminating process.\n");
                    exit(1);
                }

                uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[currentFieldIndex] - 1) - i);
                if ((fieldMask & maskBit) == 0) // wildcard bit
                {
                    if (currentNode->zeroBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->zeroBranch);
                    }

                    if (currentNode->oneBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->oneBranch);
                    }

                    if (currentNode->wildcardBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->wildcardBranch);
                    }
                } else {
                    if ((maskedFieldValue & maskBit) == 0) // zero bit
                    {
                        if (currentNode->zeroBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->zeroBranch);
                        }

                        if (currentNode->wildcardBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->wildcardBranch);
                        }
                    } else // one bit
                    {
                        if (currentNode->oneBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->oneBranch);
                        }

                        if (currentNode->wildcardBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->wildcardBranch);
                        }
                    }
                }
            }

            vCurrentLevelNodes = vNextLevelNodes;
            vNextLevelNodes.clear();
        }
        // get all next level trees
        for (auto node : vCurrentLevelNodes) {
            if (node->nextLevelTrie != nullptr) {
                vOutputTries.push_back(node->nextLevelTrie);
            } else {
                fprintf(stderr,
                        "[Trie::getNextLevelEquivalenceClasses] Invalid node (node->nextLevelTrie = nullptr) found. Terminating process.\n");
                exit(1);
            }
        }
    }

    // Found the next level tries. Now, compute the equivalence classes.

    list<uint64_t> lowerBoundList, upperBoundList;
    unordered_map<uint64_t, uint64_t> lowerBoundMap, upperBoundMap;

    fieldValue = Trie::getIntValue(nextFieldIndex, rule.fieldValue[nextFieldIndex]);
    fieldMask = Trie::getIntValue(nextFieldIndex, rule.fieldMask[nextFieldIndex]);
    maskedFieldValue = fieldValue & fieldMask;

    for (auto outputTrie : vOutputTries) {
        bool matchFound = true;
        if (outputTrie->getTotalRuleCount() == 0) {
            continue;
        }

        TrieNode *currentNode = outputTrie->root;
        EquivalenceRange range;
        for (unsigned int i = 0; i < ::fieldWidth[nextFieldIndex]; i++) {
            uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[nextFieldIndex] - 1) - i);
            if ((fieldMask & maskBit) != 0) // non-wildcard bit
            {
                if ((maskedFieldValue & maskBit) == 0) // zero bit
                {
                    if (currentNode->zeroBranch != nullptr) {
                        currentNode = currentNode->zeroBranch;
                    } else {
                        /* fprintf(stderr, "[Trie::getNextLevelEquivalenceClasses] Error: currentNode->zeroBranch cannot be nullptr. Terminating process.\n");
                        fprintf(stderr, "[Trie::getNextLevelEquivalenceClasses] currentFieldIndex: %d, lb: %lu, data: %s\n",
                                currentFieldIndex, lb, data.toString().c_str());
                        exit(1); */

                        matchFound = false;
                        break;
                    }
                } else // one bit
                {
                    range.lowerBound |= 1;
                    range.upperBound |= 1;

                    if (currentNode->oneBranch != nullptr) {
                        currentNode = currentNode->oneBranch;
                    } else {
                        /* fprintf(stderr, "[Trie::getNextLevelEquivalenceClasses] Error: currentNode->oneBranch cannot be nullptr. Terminating process.\n");
                        fprintf(stderr, "[Trie::getNextLevelEquivalenceClasses] currentFieldIndex: %d, lb: %lu, data: %s\n",
                                currentFieldIndex, lb, data.toString().c_str());
                        exit(1); */

                        matchFound = false;
                        break;
                    }
                }

                range.lowerBound <<= 1;
                range.upperBound <<= 1;
            } else // wildcard bit
            {
                matchFound = true;
                break;
            }
        }

        if (!matchFound) {
            continue;
        }

        // Perform DFS (Depth-first Search).

        stack<TrieNode *> nodes;
        nodes.push(currentNode);
        stack<EquivalenceRange> ranges;
        ranges.push(range);

        while (!nodes.empty()) {
            currentNode = nodes.top();
            nodes.pop();
            EquivalenceRange tempRange = ranges.top();
            ranges.pop();

            // fprintf(stdout, "Exploring node...\n");

            if (currentNode == nullptr) {
                fprintf(stderr,
                        "[Trie::getNextLevelEquivalenceClasses] Invalid node (node = nullptr) found. Terminating process.\n");
                exit(1);
            }

            if (currentNode->oneBranch != nullptr) {
                EquivalenceRange oneRange = tempRange;
                oneRange.lowerBound |= 1;
                oneRange.upperBound |= 1;
                oneRange.lowerBound <<= 1;
                oneRange.upperBound <<= 1;
                ranges.push(oneRange);

                // fprintf(stdout, "Pushing oneBranch...\n");
                nodes.push(currentNode->oneBranch);
            }

            if (currentNode->zeroBranch != nullptr) {
                EquivalenceRange zeroRange = tempRange;
                // zeroRange.lowerBound |= 0;
                // zeroRange.upperBound |= 0;
                zeroRange.lowerBound <<= 1;
                zeroRange.upperBound <<= 1;
                ranges.push(zeroRange);

                // fprintf(stdout, "Pushing zeroBranch...\n");
                nodes.push(currentNode->zeroBranch);
            }

            if (currentNode->wildcardBranch != nullptr) {
                EquivalenceRange wildcardRange = tempRange;
                // wildcardRange.lowerBound |= 0;
                wildcardRange.upperBound |= 1;
                wildcardRange.lowerBound <<= 1;
                wildcardRange.upperBound <<= 1;
                ranges.push(wildcardRange);

                // fprintf(stdout, "Pushing wildcardBranch...\n");
                nodes.push(currentNode->wildcardBranch);
            }

            if ((currentNode->nextLevelTrie != nullptr) || (currentNode->ruleSet != nullptr)) {
                // reached a leaf

                tempRange.lowerBound >>= 1;
                tempRange.upperBound >>= 1;

                Trie::addToBoundList(tempRange.lowerBound, tempRange.upperBound, lowerBoundList, upperBoundList,
                                     lowerBoundMap, upperBoundMap);
            } else {
                // Do nothing.
            }
        }
    }

    lowerBoundList.sort();
    upperBoundList.sort();

    // fprintf(stdout, "Generating equivalent packets...\n");

    uint64_t lowerBound = 0;
    uint64_t upperBound = 0;

    if (!lowerBoundList.empty()) {
        lowerBound = upperBound = lowerBoundList.front();
    }

    while (!upperBoundList.empty()) {
        if (!lowerBoundList.empty()) {
            lowerBound = lowerBoundList.front();
            if ((lowerBound > (upperBound + 1))
                /* && (upperBoundList.front() > (upperBound + 1))
                && (upperBoundMap[upperBoundList.front()] < (upperBound + 1)) */) {
                lowerBound = upperBound + 1;

                if (lowerBoundList.front() <= upperBoundList.front()) {
                    upperBound = lowerBoundList.front() - 1;
                } else {
                    upperBound = upperBoundList.front();
                    upperBoundList.pop_front();
                }

                EquivalenceClass packetClass;
                packetClass.lowerBound[nextFieldIndex] = lowerBound;
                packetClass.upperBound[nextFieldIndex] = upperBound;

                vPacketClasses.push_back(packetClass);

                continue;
            } else {
                lowerBoundList.pop_front();
            }
        } else {
            lowerBound = upperBound + 1;
        }

        if (!lowerBoundList.empty()) {
            if (lowerBoundList.front() <= upperBoundList.front()) {
                upperBound = lowerBoundList.front() - 1;
            } else {
                upperBound = upperBoundList.front();
                upperBoundList.pop_front();
            }
        } else {
            upperBound = upperBoundList.front();
            upperBoundList.pop_front();
        }

        EquivalenceClass packetClass;
        packetClass.lowerBound[nextFieldIndex] = lowerBound;
        packetClass.upperBound[nextFieldIndex] = upperBound;

        vPacketClasses.push_back(packetClass);
    }
}

/*
 * Only the last level trie is used for each equivalence class to build the corresponding forwarding graph.
 * Ignore the DUMMY rules, or only put the FORWARDING rules in the graph.
 * */
ForwardingGraph *Trie::getForwardingGraph(FieldIndex currentFieldIndex, const vector<Trie *> &vInputTries,
                                          const EquivalenceClass &packetClass, FILE *fp) {
    ForwardingGraph *graph = new ForwardingGraph;
    if (graph == nullptr) {
        fprintf(stderr, "[Trie::getForwardingGraph] Memory allocation error (graph == nullptr). Terminating process.\n");
        exit(1);
    }

    if ((currentFieldIndex + 1) != ALL_FIELD_INDEX_END_MARKER) {
        fprintf(stderr,
                "[Trie::getForwardingGraph] Method called on wrong trie (field index = %d). Terminating process.\n",
                currentFieldIndex);
        exit(1);
    }

    uint64_t fieldValue = packetClass.lowerBound[currentFieldIndex];
    uint64_t fieldMask = EquivalenceClass::getMaxValue(currentFieldIndex);

    uint64_t maskedFieldValue = fieldValue & fieldMask;

    for (unsigned int t = 0; t < vInputTries.size(); t++) {
        Trie *inputTrie = vInputTries[t];
        if (inputTrie->getTotalRuleCount() == 0) {
            continue;
        }

        vector<TrieNode *> vCurrentLevelNodes, vNextLevelNodes;
        vCurrentLevelNodes.push_back(inputTrie->root);
        TrieNode *currentNode = nullptr;

        for (unsigned int i = 0; i < ::fieldWidth[currentFieldIndex]; i++) {
            while (!vCurrentLevelNodes.empty()) {
                currentNode = vCurrentLevelNodes.at(0);
                vCurrentLevelNodes.erase(vCurrentLevelNodes.begin());

                if (currentNode == nullptr) {
                    fprintf(stderr,
                            "[Trie::getForwardingGraph] Invalid node (node = nullptr) found. Terminating process.\n");
                    exit(1);
                }

                uint64_t maskBit = (uint64_t) 1 << ((unsigned int) (::fieldWidth[currentFieldIndex] - 1) - i);
                if ((fieldMask & maskBit) == 0) // wildcard bit
                {
                    if (currentNode->zeroBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->zeroBranch);
                    }

                    if (currentNode->oneBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->oneBranch);
                    }

                    if (currentNode->wildcardBranch != nullptr) {
                        vNextLevelNodes.push_back(currentNode->wildcardBranch);
                    }
                } else {
                    if ((maskedFieldValue & maskBit) == 0) // zero bit
                    {
                        if (currentNode->zeroBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->zeroBranch);
                        }

                        if (currentNode->wildcardBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->wildcardBranch);
                        }
                    } else // one bit
                    {
                        if (currentNode->oneBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->oneBranch);
                        }

                        if (currentNode->wildcardBranch != nullptr) {
                            vNextLevelNodes.push_back(currentNode->wildcardBranch);
                        }
                    }
                }
            }

            vCurrentLevelNodes = vNextLevelNodes;
            vNextLevelNodes.erase(vNextLevelNodes.begin(), vNextLevelNodes.end());
        }

        for (unsigned int i = 0; i < vCurrentLevelNodes.size(); i++) {
            TrieNode *node = vCurrentLevelNodes[i];
            if (node->ruleSet != nullptr) {
                unordered_set<Rule, KHash<Rule>, KEqual<Rule> >::const_iterator itr;
                for (itr = node->ruleSet->begin(); itr != node->ruleSet->end(); ++itr) {
                    const Rule &rule = *itr;
                    if (rule.type != FORWARDING) {
                        continue;
                    }

                    if (rule.priority == INVALID_PRIORITY) {
                        continue;
                    }

                    ForwardingLink link(rule, false);

                    if (mode == TEST_MODE) {
                        // For the testVerification() experiment present in Test.cpp.
                        if (rule.location.compare(rule.nextHop) == 0) {
                            link.isGateway = true;
                        }
                    } else if (mode == PROXY_MODE) {
                        if (rule.nextHop.compare(rule.fieldValue[NW_DST]) == 0) {
                            link.isGateway = true;
                        }
                        for (unsigned int i = 0; i < endhosts.size(); i++) {
                            if (rule.nextHop.compare(endhosts[i]) == 0) {
                                link.isGateway = true;
                                break;
                            }
                        }
                    }

                    // fprintf(stdout, "[Trie::getForwardingGraph] %s\n", link.toString().c_str());
                    // fprintf(fp, "[Trie::getForwardingGraph] %s\n", link.toString().c_str());
                    graph->addLink(link);
                }
            } else {
                fprintf(stderr,
                        "[Trie::getForwardingGraph] Invalid node (node->ruleSet = nullptr) found. Terminating process.\n");
                exit(1);
            }
        }
    }

    return graph;
}

int Trie::getTotalRuleCount() const {
    return this->totalRuleCount;
}

void Trie::print(FILE *fp) const {
    Trie::traversePreorder(this->root, this->fieldIndex, 0, -1, fp);
}

void Trie::traversePreorder(const TrieNode *node, FieldIndex index, int level, int branch, FILE *fp) {
    if (node == nullptr) {
        return;
    }

    string linePrefix = "";
    string leadingSpace = "";

    for (int i = 0; i < level; i++) {
        linePrefix += "-";
        leadingSpace += " ";
    }

    fprintf(fp, "%s", linePrefix.c_str());

    if (branch == -1) {
        fprintf(fp, "[root] Field index: %d", index);
    } else {
        fprintf(fp, "%d", branch);
    }

    if (node->ruleSet != nullptr) {
        unordered_set<Rule, KHash<Rule>, KEqual<Rule> >::const_iterator itr;
        for (itr = node->ruleSet->begin(); itr != node->ruleSet->end(); ++itr) {
            fprintf(fp, "\n");
            fprintf(fp, "%s", leadingSpace.c_str());
            fprintf(fp, "%s", itr->toString().c_str());
        }
    }

    fprintf(fp, "\n");

    if (node->nextLevelTrie != nullptr) {
        node->nextLevelTrie->print(fp);
    }

    Trie::traversePreorder(node->zeroBranch, index, level + 1, 0, fp);
    Trie::traversePreorder(node->oneBranch, index, level + 1, 1, fp);
    Trie::traversePreorder(node->wildcardBranch, index, level + 1, 2, fp);
}
