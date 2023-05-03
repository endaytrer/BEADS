/*
 * VeriFlow.cpp
 *
 *  Created on: Mar 12, 2012
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
#include <fstream>
#include <string>
#include <csignal>
#include "net.h"
#include "thread.h"
#include "StringTokenizer.h"
#include "VeriFlow.h"
#include "openflow.h"
#include "OpenFlowProtocolMessage.h"
#include "DynamicArray.h"
#include "Network.h"
#include <list>
#include <vector>
#include <sys/time.h>
#include <unordered_map>
#include <unordered_set>
#include "EquivalenceClass.h"
#include "Rule.h"
#include "ForwardingGraph.h"
#include "ForwardingLink.h"
#include "Template.h"
#include "Trie.h"
#include "Test.h"

using namespace std;

static int tcpServerSocket;

static string controllerIpAddress = "127.0.0.1"; // localhost
static unsigned short controllerPort = 6633; // default NOX port

static FILE *logFile = nullptr;

static pthread_mutex_t networkMutex, veriflowMutex;

static VeriFlow veriflow;
Network network;

int mode = TEST_MODE;

vector<string> endhosts;

vector<EquivalenceClass> faults;

int main(int argc, char **argv) {
    if (argc == 1) {
        mode = TEST_MODE;
        Test::test();
        return EXIT_SUCCESS;
    } else if (argc < 6) {
        fprintf(stderr, "USAGE: %s <veriflow_port> <controller_address> <controller_port> <topology_file> <log_file>\n",
                argv[0]);
        exit(1);
    }

    mode = PROXY_MODE;

    // Network network;

    string topologyFileName = argv[4];
    parseTopologyFile(topologyFileName, network);

    network.print();
    // return EXIT_SUCCESS;

    controllerIpAddress = argv[2];
    controllerPort = (unsigned short) atoi(argv[3]);

    tcpServerSocket = createSocket(SOCK_STREAM);
    logFile = fopen(argv[5], "w");

    struct sigaction act{};
    act.sa_handler = signalHandler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = 0;
    sigaction(SIGINT, &act, nullptr);

    int optval = 1;

    setsockopt(tcpServerSocket, SOL_SOCKET, SO_REUSEADDR, &optval, sizeof(optval));

    bindSocket(tcpServerSocket, nullptr, (unsigned short) atoi(argv[1]));
    listenSocket(tcpServerSocket, BACKLOG);

    createMutex(&networkMutex);
    createMutex(&veriflowMutex);

    unsigned int i = 0;

    while (true) {
        struct sockaddr_in clientAddress{};
        socklen_t clientAddressLength = sizeof(clientAddress);
        int clientSocket;

        ++i;

        fprintf(stdout, "[VeriFlow] [%u] Waiting for new connection...\n", i);
        clientSocket = accept(tcpServerSocket, (struct sockaddr *) &clientAddress, &clientAddressLength);
        fprintf(stdout, "[VeriFlow] [%u] Accepted new connection from %s at port %u\n",
                i, getIPAddress(clientAddress), ntohs(clientAddress.sin_port));

        VeriFlowConnectionInfo info{};
        info.clientSocket = clientSocket;
        info.clientAddress = clientAddress;
        info.network = &network;
        info.veriflow = &veriflow;
        info.networkMutex = &networkMutex;
        info.veriflowMutex = &veriflowMutex;

        handleVeriFlowConnection(info);
    }

    return EXIT_SUCCESS;
}

void parseTopologyFile(const string &fileName, Network &io_network) {
    char buffer[1024];
    ifstream fin(fileName.c_str());
    while (!fin.eof()) {
        fin.getline(buffer, 1023);
        if (strlen(buffer) == 0) {
            continue;
        }

        if (strstr(buffer, "#") == buffer) {
            continue;
        }

        // # Format: id ipAddress endDevice port1 nextHopIpAddress1 port2 nextHopIpAddress2 ...

        StringTokenizer st(buffer, " ");
        if (st.countTokens() < 3) {
            continue;
        }

        unsigned int id = atoi(st.nextToken().c_str());
        string ipAddress = st.nextToken();
        bool endDevice = atoi(st.nextToken().c_str());

        if (endDevice) {
            endhosts.push_back(ipAddress);
        }

        io_network.addDevice(id, ipAddress, endDevice);

        while (st.hasMoreTokens()) {
            string port = st.nextToken();
            if (port == "#") {
                break;
            }

            string nextHopIpAddress = st.nextToken();

            io_network.addPort(ipAddress, atoi(port.c_str()), nextHopIpAddress);
        }
    }

    fin.close();
}

void handleVeriFlowConnection(VeriFlowConnectionInfo &info) {
    // Connect to the controller.
    int controllerSocket = createSocket(SOCK_STREAM);
    struct sockaddr_in controllerAddress = createSocketAddress(controllerIpAddress.c_str(), controllerPort);

    int res = connect(controllerSocket, (struct sockaddr *) &controllerAddress, sizeof(controllerAddress));
    if (res == -1) {
        fprintf(stderr, "[handleVeriFlowConnection] Cannot connect to controller at %s.\n",
                controllerIpAddress.c_str());
        close(controllerSocket);
        close(info.clientSocket);

        return;
    }

    info.controllerSocket = controllerSocket;

    ProxyConnectionInfo *info1 = new ProxyConnectionInfo;
    info1->clientAddress = info.clientAddress;
    info1->recvSocket = info.controllerSocket;
    info1->sendSocket = info.clientSocket;
    info1->network = info.network;
    info1->veriflow = info.veriflow;
    info1->networkMutex = info.networkMutex;
    info1->veriflowMutex = info.veriflowMutex;

    ProxyConnectionInfo *info2 = new ProxyConnectionInfo;
    info2->clientAddress = info.clientAddress;
    info2->recvSocket = info.clientSocket;
    info2->sendSocket = info.controllerSocket;
    info2->network = info.network;
    info2->veriflow = info.veriflow;
    info2->networkMutex = info.networkMutex;
    info2->veriflowMutex = info.veriflowMutex;

    info1->other = info2;
    info2->other = info1;

    pthread_t controllerToNetworkCommunicationThread;
    createThread(&controllerToNetworkCommunicationThread, proxyCommunicationThreadFunction, (void *) info1,
                 PTHREAD_CREATE_DETACHED, NORMAL_PRIORITY);

    pthread_t networkToControllerCommunicationThread;
    createThread(&networkToControllerCommunicationThread, proxyCommunicationThreadFunction, (void *) info2,
                 PTHREAD_CREATE_DETACHED, NORMAL_PRIORITY);
}

void *proxyCommunicationThreadFunction(void *arg) {
    setThreadAsyncCancel();

    ProxyConnectionInfo *info = (ProxyConnectionInfo *) arg;

    DynamicArray<char> messageBuffer;
    char data[MAX_BUFFER_SIZE];
    int bytesReceived = 0;

    while ((bytesReceived = recv(info->recvSocket, data, sizeof(data), 0)) > 0) {
/*		int bytesSent = 0, res = 0;

		while(bytesSent < bytesReceived)
		{
			res = send(info.sendSocket, (data + bytesSent), (bytesReceived - bytesSent), 0);
			if(res == -1)
			{
				fprintf(stderr, "[proxyCommunicationThreadFunction] TCP send failure. Stopping operation.\n");

				close(info.sendSocket);
				close(info.recvSocket);

				pthread_exit(nullptr);
			}
			else
			{
				bytesSent += res;
			}
		}

		continue;
*/
        messageBuffer.append(data, bytesReceived);

        while (messageBuffer.size() >= sizeof(ofp_header)) {
            // Received full header.
            const ofp_header *header = (const ofp_header *) messageBuffer.getData(sizeof(ofp_header));

            if (messageBuffer.size() >= ntohs(header->length)) {
                // Received a complete message.
                int bytesToSend = ntohs(header->length);
                char *messageData = messageBuffer.getData(bytesToSend);

                delete[] (char *) header;

                // pthread_mutex_lock(info.mutex);
                // fprintf(logFile, "[%s]\n", getIPAddress(info.clientAddress));
                OpenFlowProtocolMessage::process(messageData, *info, logFile);
                // fprintf(logFile, "\n");
                // fflush(logFile);
                // pthread_mutex_unlock(info.mutex);

                int bytesSent = 0, res = 0;

                while (bytesSent < bytesToSend) {
                    res = send(info->sendSocket, (messageData + bytesSent), (bytesToSend - bytesSent), 0);
                    if (res == -1) {
                        fprintf(stderr, "[proxyCommunicationThreadFunction] TCP send failure. Stopping operation.\n");

                        info->other->recvSocket = -1;
                        close(info->sendSocket);
                        info->other->sendSocket = -1;
                        close(info->recvSocket);

                        pthread_exit(nullptr);
                    } else {
                        bytesSent += res;
                    }
                }

                delete[] messageData;

                messageBuffer.clearRange(0, bytesToSend - 1);
            } else {
                break;
            }
        }
    }

    if (bytesReceived == 0) {
        info->other->sendSocket = -1;
        close(info->recvSocket);
        info->other->recvSocket = -1;
        close(info->sendSocket);
    }

    fprintf(stdout, "[proxyCommunicationThreadFunction] Connection closed.\n");

    pthread_exit(nullptr);
}

void signalHandler(int sig) {
    close(tcpServerSocket);

    pthread_mutex_lock(&veriflowMutex);
    //veriflow.print(logFile);
    fclose(logFile);
    pthread_mutex_unlock(&veriflowMutex);

    exit(1);
}

uint64_t getMacValueAsInt(const string &macAddress) {
    uint64_t macValue = 0;
    StringTokenizer st(macAddress, ":");
    while (st.hasMoreTokens()) {
        unsigned long value = strtoul(st.nextToken().c_str(), nullptr, 16);
        macValue <<= 8;
        macValue += value;
    }

    return macValue;
}

string getMacValueAsString(uint64_t macAddress) {
    unsigned int values[6];

    for (int i = 5; i >= 0; i--) {
        values[5 - i] = (unsigned int) ((macAddress >> (8 * i)) & ((unsigned int) 0xFF));
    }

    char buffer[1024];
    sprintf(buffer, "%02x:%02x:%02x:%02x:%02x:%02x", values[0], values[1], values[2], values[3], values[4], values[5]);

    string macValue = buffer;
    return macValue;
}

string getMacValueAsString(const uint8_t *macAddress) {
    string macValue;
    char buffer[8];

    for (int i = 0; i < OFP_ETH_ALEN; i++) {
        uint8_t upperNibble = macAddress[i] >> 4;
        uint8_t lowerNibble = macAddress[i] & (uint8_t) 0xF;

        sprintf(buffer, "%x%x", upperNibble, lowerNibble);
        macValue += buffer;
        if (i != (OFP_ETH_ALEN - 1)) {
            macValue += ":";
        }
    }

    return macValue;
}

uint64_t getIpValueAsInt(const string &ipAddress) {
    uint64_t ipValue = 0;
    StringTokenizer st(ipAddress, ".");
    while (st.hasMoreTokens()) {
        unsigned int quadValue = atoi(st.nextToken().c_str());
        ipValue <<= 8;
        ipValue += quadValue;
    }

    return ipValue;
}

string getIpValueAsString(uint64_t ipAddress) {
    unsigned int quadValues[4];

    for (int i = 3; i >= 0; i--) {
        quadValues[3 - i] = (unsigned int) ((ipAddress >> (8 * i)) & ((unsigned int) 0xFF));
    }

    char buffer[1024];
    sprintf(buffer, "%u.%u.%u.%u", quadValues[0], quadValues[1], quadValues[2], quadValues[3]);

    string ipValue = buffer;
    return ipValue;
}

string convertMaskToDottedFormat(unsigned int mask) {
    uint64_t maskValue = 0xFFFFFFFF;
    maskValue >>= mask;
    maskValue <<= mask;

    return ::getIpValueAsString(maskValue);
}

string convertIntToString(unsigned int value) {
    char buffer[128];
    snprintf(buffer, 127, "%u", value);

    string retVal = buffer;
    return retVal;
}

bool compareForwardingLink(const ForwardingLink &first, const ForwardingLink &second) {
    if (first.rule.priority >= second.rule.priority) {
        return true;
    } else {
        return false;
    }
}

VeriFlow::VeriFlow() {
    this->primaryTrie = new Trie(IN_PORT);
    this->previousFailures = 0;
    if (this->primaryTrie == nullptr) {
        fprintf(stderr,
                "[VeriFlow::VeriFlow] Memory allocation error (this->primaryTrie == nullptr). Terminating process.\n");
        exit(1);
    }
}

VeriFlow::~VeriFlow() {
    if (this->primaryTrie != nullptr) {
        delete this->primaryTrie;
        this->primaryTrie = nullptr;
    }
}

bool VeriFlow::addRule(const Rule &rule) {
    // add rule to every layers in the tree

    Trie *currentTrie = this->primaryTrie;
    vector<Trie *> vTries;
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        vTries.push_back(currentTrie);
        TrieNode *leaf = currentTrie->addRule(rule);
        if (i == (ALL_FIELD_INDEX_END_MARKER - 1)) {
            // This was the last level trie. Need to check the rule list.
            if (leaf->ruleSet == nullptr) {
                leaf->ruleSet = new unordered_set<Rule, KHash<Rule>, KEqual<Rule> >;
                if (leaf->ruleSet == nullptr) {
                    fprintf(stderr,
                            "[VeriFlow::addRule] Memory allocation error (leaf->ruleSet == nullptr). Terminating process.\n");
                    exit(1);
                }
            } else {
                unordered_set<Rule, KHash<Rule>, KEqual<Rule> >::const_iterator itr;
                itr = leaf->ruleSet->find(rule);
                if (itr != leaf->ruleSet->end()) // Rule already exists.
                {
                    return false;
                }
            }

            leaf->ruleSet->insert(rule);
        } else {
            // This was an intermediate trie. Only add rules to the last level trees.
            if (leaf->nextLevelTrie == nullptr) {
                leaf->nextLevelTrie = new Trie((FieldIndex) (i + 1));
                if (leaf->nextLevelTrie == nullptr) {
                    fprintf(stderr,
                            "[VeriFlow::addRule] Memory allocation error (leaf->nextLevelTrie == nullptr). Terminating process.\n");
                    exit(1);
                }
            }

            currentTrie = leaf->nextLevelTrie;
        }
    }

    for (auto & vTrie : vTries) {
        (vTrie->totalRuleCount)++;
    }

    return true;
}

bool VeriFlow::removeRule(const Rule &rule) {
    Trie *currentTrie = this->primaryTrie;
    vector<Trie *> vTries;
    vector<TrieNode *> vLeaves;
    for (int i = 0; i < ALL_FIELD_INDEX_END_MARKER; i++) {
        TrieNode *leaf = currentTrie->findNode(rule.fieldValue[i], rule.fieldMask[i]);
        if (leaf == nullptr) {
            return false;
        }

        if (i == (ALL_FIELD_INDEX_END_MARKER - 1)) {
            // This was the last level trie. Need to check the rule list.
            if (leaf->ruleSet == nullptr) {
                fprintf(stderr, "[VeriFlow::removeRule] Error: leaf->ruleSet cannot be nullptr. Terminating process.\n");
                exit(1);
            }

            unordered_set<Rule, KHash<Rule>, KEqual<Rule> >::const_iterator itr;
            itr = leaf->ruleSet->find(rule);
            if (itr != leaf->ruleSet->end()) // Rule found.
            {
                leaf->ruleSet->erase(itr);
                if (leaf->ruleSet->empty()) {
                    currentTrie->removeRule(leaf);

                    for (unsigned int k = 0; k < vLeaves.size(); k++) {
                        unsigned int index = (vLeaves.size() - k - 1);
                        if (vLeaves[index]->nextLevelTrie->totalRuleCount == 0) {
                            vTries[index]->removeRule(vLeaves[index]);
                        }
                    }
                }

                return true;
            }

            return false;
        } else {
            // This was an intermediate trie.

            vTries.push_back(currentTrie);
            vLeaves.push_back(leaf);

            if (leaf->nextLevelTrie == nullptr) {
                fprintf(stderr,
                        "[VeriFlow::removeRule] Error: leaf->nextLevelTrie cannot be nullptr. Terminating process.\n");
                exit(1);
            }

            currentTrie = leaf->nextLevelTrie;
        }
    }

    return false;
}

void
VeriFlow::traverseAllTries(const Rule &rule, const vector<Trie *> &prevTries, const vector<EquivalenceClass> &prevClasses, FieldIndex fi, uint64_t *lb, uint64_t *ub, vector<EquivalenceClass> &vFinalPacketClasses, vector<vector<Trie *> > &vFinalTries) {
    if (fi == ALL_FIELD_INDEX_END_MARKER - 1) {
        for (auto &prevClass: prevClasses) {
            lb[fi] = prevClass.lowerBound[fi];
            ub[fi] = prevClass.upperBound[fi];
            EquivalenceClass packetClass(
                    lb,
                    ub);

            vFinalPacketClasses.push_back(
                    packetClass);
            vFinalTries.push_back(
                    prevTries);
        }
        return;
    }
    for (auto &prevClass: prevClasses) {
        if (rule.type != FORWARDING) {
            continue;
        }

        lb[fi] = prevClass.lowerBound[fi];
        ub[fi] = prevClass.upperBound[fi];
        vector<EquivalenceClass> nextClasses;
        vector<Trie *> nextTries;
        Trie::getNextLevelEquivalenceClasses(fi, prevClass.lowerBound[fi], rule, prevTries, nextClasses, nextTries);
        if (nextClasses.empty()) {
            fprintf(stderr, "[VeriFlow::getAffectedEquivalenceClasses] Error in rule: %s\n",
                    rule.toString().c_str());
            fprintf(stderr,
                    "[VeriFlow::getAffectedEquivalenceClasses] Error: ((FieldIndex #%d's equivalence class).size() == 0). Terminating process.\n", fi + 1);
            exit(1);
        }
        VeriFlow::traverseAllTries(rule, nextTries, nextClasses, (FieldIndex)((int)fi + 1), lb, ub, vFinalPacketClasses, vFinalTries);
    }
}

bool
VeriFlow::getAffectedEquivalenceClasses(const Rule &rule, int command, vector<EquivalenceClass> &vFinalPacketClasses,
                                        vector<vector<Trie *> > &vFinalTries) {
    if (command == OFPFC_ADD) {
        // We may choose not to verify a rule if that rule is already present in the data plane.
        bool res = this->addRule(rule);
        if (!res) {
            return false;
        }
    } else if (command == OFPFC_DELETE_STRICT) {
        bool res = this->removeRule(rule);
        if (!res) {
            return false;
        } else {
            Rule dummyRule = rule;
            dummyRule.type = DUMMY;
            this->addRule(dummyRule); // This dummy rule will be deleted inside VeriFlow::verifyRule() function.
        }
    }

    vector<EquivalenceClass> vInPortPacketClasses;
    this->primaryTrie->getEquivalenceClasses(rule, vInPortPacketClasses);

    if (vInPortPacketClasses.empty()) {
        fprintf(stderr, "[VeriFlow::getAffectedEquivalenceClasses] Error in rule: %s\n", rule.toString().c_str());
        fprintf(stderr,
                "[VeriFlow::getAffectedEquivalenceClasses] Error: (vInPortPacketClasses.size() == 0). Terminating process.\n");
        exit(1);
    }

    vector<Trie *> vInPortTries;
    vInPortTries.push_back(this->primaryTrie);
    uint64_t lb[ALL_FIELD_INDEX_END_MARKER], ub[ALL_FIELD_INDEX_END_MARKER];
    traverseAllTries(rule, vInPortTries, vInPortPacketClasses, IN_PORT, lb, ub, vFinalPacketClasses, vFinalTries);
    return true;
}

void
VeriFlow::processCurrentHop(const EquivalenceClass &packetClass, ForwardingGraph *graph, const string &currentLocation,
                            unordered_set<string> &visited, NextHopInfo &nextHopInfo, FILE *fp) {
    if (graph == nullptr) {
        // fprintf(fp, "[VeriFlow::processCurrentHop] (graph == nullptr) for the following packet class.\n");
        // fprintf(fp, "[VeriFlow::processCurrentHop] PacketClass: %s\n", packetClass.toString().c_str());
        return;
    }

    if (visited.find(currentLocation) != visited.end()) {
        // Found a loop.
        // fprintf(fp, "[VeriFlow::processCurrentHop] Found a loop for the following packet class.\n");
        // fprintf(fp, "[VeriFlow::processCurrentHop] PacketClass: %s\n", packetClass.toString().c_str());

        return;
    }

    if (graph->links.find(currentLocation) == graph->links.end()) {
        // Found a black hole.
        // fprintf(fp, "[VeriFlow::processCurrentHop] Found a black hole for the following packet class as current location (%s) not found.\n", currentLocation.c_str());
        // fprintf(fp, "[VeriFlow::processCurrentHop] PacketClass: %s\n", packetClass.toString().c_str());

        return;
    }

    if (graph->links[currentLocation].empty()) {
        // Found a black hole.
        // fprintf(fp, "[VeriFlow::processCurrentHop] Found a black hole for the following packet class as there is no outgoing link at current location (%s).\n", currentLocation.c_str());
        // fprintf(fp, "[VeriFlow::processCurrentHop] PacketClass: %s\n", packetClass.toString().c_str());

        return;
    }

    graph->links[currentLocation].sort(compareForwardingLink);

    const list<ForwardingLink> &linkList = graph->links[currentLocation];
    list<ForwardingLink>::const_iterator itr = linkList.begin();

    nextHopInfo.nextHop = itr->rule.nextHop;
    nextHopInfo.visited = visited;
    nextHopInfo.visited.insert(currentLocation);
}

bool VeriFlow::verifyRule(const Rule &rule, int command, double &updateTime, double &packetClassSearchTime,
                          double &graphBuildTime, double &queryTime, unsigned long &ecCount, FILE *fp) {
    // fprintf(fp, "[VeriFlow::verifyRule] verifying this rule: %s\n", rule.toString().c_str());

    updateTime = packetClassSearchTime = graphBuildTime = queryTime = 0;
    ecCount = 0;

    struct timeval start, end;
    double usecTime, seconds, useconds;

    gettimeofday(&start, nullptr);
    // May add code in a future version to maintain a cache of forwarding graphs. This cache needs to be updated for every new rule.
    gettimeofday(&end, nullptr);

    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    usecTime = (seconds * 1000000) + useconds;
    updateTime = usecTime;

    gettimeofday(&start, nullptr);
    vector<EquivalenceClass> vFinalPacketClasses;
    vector<vector<Trie *> > vFinalTries;

    bool res = this->getAffectedEquivalenceClasses(rule, command, vFinalPacketClasses, vFinalTries);
    if (!res) {
        return false;
    }
    gettimeofday(&end, nullptr);

    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    usecTime = (seconds * 1000000) + useconds;
    packetClassSearchTime = usecTime; // count time

    ecCount = vFinalPacketClasses.size();
    if (ecCount == 0) {
        fprintf(stderr, "[VeriFlow::verifyRule] Error in rule: %s\n", rule.toString().c_str());
        fprintf(stderr,
                "[VeriFlow::verifyRule] Error: (ecCount = vFinalPacketClasses.size() = 0). Terminating process.\n");
        exit(1);
    } else {
        // fprintf(stdout, "\n");
        // fprintf(stdout, "[VeriFlow::verifyRule] ecCount: %lu\n", ecCount);
    }

    // fprintf(stdout, "[VeriFlow::verifyRule] Generating forwarding graphs...\n");
    gettimeofday(&start, nullptr);
    vector<ForwardingGraph *> vGraph;
    for (unsigned int i = 0; i < vFinalPacketClasses.size(); i++) {
        EquivalenceClass packetClass = vFinalPacketClasses[i];
        // fprintf(stdout, "[VeriFlow::verifyRule] [%u] ecCount: %lu, %s\n", i, ecCount, packetClass.toString().c_str());
        ForwardingGraph *graph = Trie::getForwardingGraph(TP_DST, vFinalTries[i], packetClass, fp);
        vGraph.push_back(graph);
    }
    gettimeofday(&end, nullptr);
    // fprintf(stdout, "[VeriFlow::verifyRule] Generated forwarding graphs.\n");

    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    usecTime = (seconds * 1000000) + useconds;
    graphBuildTime = usecTime;

    // fprintf(stdout, "[VeriFlow::verifyRule] Running query...\n");
    gettimeofday(&start, nullptr);
    // Add query code here
    size_t currentFailures = 0;
    for (unsigned int i = 0; i < vGraph.size(); i++) {
        unordered_set<string> visited;
        string lastHop = network.getNextHopIpAddress(rule.location, rule.in_port);
        // From one step before the designated string
        // fprintf(fp, "start traversing at: %s\n", rule.location.c_str());
        if (!this->traverseForwardingGraph(vFinalPacketClasses[i], vGraph[i], rule.location, lastHop, visited, fp)) {
            ++currentFailures;
        }
    }

    fprintf(stderr, "faults size: %zu\n", faults.size());
    if (previousFailures > 0 && faults.empty()) {
        fprintf(fp, "[Veriflow::verifyRule] Network Fixed!\n");
    } else if (previousFailures == 0 && faults.size() > 0) {
        fprintf(fp, "[Veriflow::verifyRule] Network Broken!\n");
    }
    fflush(fp);
    previousFailures = faults.size();
    // fprintf(stdout, "[VeriFlow::verifyRule] Query complete.\n");

    if (command == OFPFC_ADD) {
        // Do nothing.
    } else if (command == OFPFC_DELETE_STRICT) {
        Rule dummyRule = rule;
        dummyRule.type = DUMMY;

        this->removeRule(dummyRule);
    }
    gettimeofday(&end, nullptr);

    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    usecTime = (seconds * 1000000) + useconds;
    queryTime = usecTime;

    for (unsigned int i = 0; i < vGraph.size(); i++) {
        delete vGraph[i];
    }

    return true;
}

bool VeriFlow::traverseForwardingGraph(const EquivalenceClass &packetClass, ForwardingGraph *graph,
                                       const string &currentLocation, const string &lastHop,
                                       unordered_set<string> visited, FILE *fp) {

    // fprintf(fp, "traversing at node: %s\n", currentLocation.c_str());
    if (graph == nullptr) {
        /* fprintf(fp, "\n");
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] (graph == nullptr) for the following packet class at node %s.\n", currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str()); */

        return true;
    }

    if (currentLocation.empty()) {
        return true;
    }

    if (visited.find(currentLocation) != visited.end()) {
        // Node is visited, Found a loop.
        fprintf(fp, "\n");
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] Found a LOOP for the following packet class at node %s.\n",
                currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str());
        for (unsigned int i = 0; i < faults.size(); i++) {
            if (packetClass.subsumes(faults[i])) {
                faults.erase(faults.begin() + i);
                i--;
            }
        }
        faults.push_back(packetClass);

        return false;
    }

    visited.insert(currentLocation);

    if (graph->links.find(currentLocation) == graph->links.end()) {
        // Found a black hole.
        fprintf(fp, "\n");
        fprintf(fp,
                "[VeriFlow::traverseForwardingGraph] Found a BLACK HOLE for the following packet class as current location (%s) not found in the graph.\n",
                currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str());
        for (unsigned int i = 0; i < faults.size(); i++) {
            if (packetClass.subsumes(faults[i])) {
                faults.erase(faults.begin() + i);
                i--;
            }
        }
        faults.push_back(packetClass);

        return false;
    }

    if (graph->links[currentLocation].empty()) {
        // Found a black hole.
        fprintf(fp, "\n");
        fprintf(fp,
                "[VeriFlow::traverseForwardingGraph] Found a BLACK HOLE for the following packet class as there is no outgoing link at current location (%s).\n",
                currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str());
        for (unsigned int i = 0; i < faults.size(); i++) {
            if (packetClass.subsumes(faults[i])) {
                faults.erase(faults.begin() + i);
                i--;
            }
        }
        faults.push_back(packetClass);

        return false;
    }

    graph->links[currentLocation].sort(compareForwardingLink); // sort using priority, and forward

    const list<ForwardingLink> &linkList = graph->links[currentLocation];
    list<ForwardingLink>::const_iterator itr = linkList.begin();
    // input_port as a filter
    if (lastHop == "nullptr" || itr->rule.in_port == 65536) {
        // do nothing
    } else {
        while (itr != linkList.end()) {
            string connected_hop = network.getNextHopIpAddress(currentLocation, itr->rule.in_port);
            if (connected_hop == lastHop) break;
            itr++;
        }
    }

    if (itr == linkList.end()) {
        // Found a black hole.
        //QUESTION: Why this is a blackhole?
        fprintf(fp, "\n");
        fprintf(fp,
                "[VeriFlow::traverseForwardingGraph] Found a BLACK HOLE for the following packet class as there is no outgoing link at current location (%s).\n",
                currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str());
        for (unsigned int i = 0; i < faults.size(); i++) {
            if (packetClass.subsumes(faults[i])) {
                faults.erase(faults.begin() + i);
                i--;
            }
        }
        faults.push_back(packetClass);

        return false;
    }

    if (itr->isGateway) {
        // Destination reachable.
        // fprintf(fp, "[VeriFlow::traverseForwardingGraph] Destination reachable.\n");
        fprintf(fp, "\n");
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] The following packet class reached destination at node %s.\n",
                currentLocation.c_str());
        fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str());
        for (unsigned int i = 0; i < faults.size(); i++) {
            // If previous fault is subset of the rule, remove the rule
            if (packetClass.subsumes(faults[i])) {
                fprintf(stderr, "Removing fault!\n");
                faults.erase(faults.begin() + i);
                i--;
            }
        }
        return true;
    } else {
        // Move to the next location.
        // fprintf(fp, "[VeriFlow::traverseForwardingGraph] Moving to node %s.\n", itr->rule.nextHop.c_str());

        if (itr->rule.nextHop == "") {
            // This rule is a packet filter. It drops packets.
            /* fprintf(fp, "\n");
            fprintf(fp, "[VeriFlow::traverseForwardingGraph] The following packet class is dropped by a packet filter at node %s.\n", currentLocation.c_str());
            fprintf(fp, "[VeriFlow::traverseForwardingGraph] PacketClass: %s\n", packetClass.toString().c_str()); */
        }

        return this->traverseForwardingGraph(packetClass, graph, itr->rule.nextHop, currentLocation, visited, fp);
    }
}

int VeriFlow::getTotalRuleCount() const {
    return this->primaryTrie->getTotalRuleCount();
}

void VeriFlow::print(FILE *fp) const {
    this->primaryTrie->print(fp);
}

void VeriFlow::setDatapathId(unsigned short socketPort, uint64_t datapathId) {
    this->socketPortToDatapathIdMap[socketPort] = datapathId;
}

uint64_t VeriFlow::getDatapathId(unsigned short socketPort) {
    if (this->socketPortToDatapathIdMap.find(socketPort) == this->socketPortToDatapathIdMap.end()) {
        return 0;
    } else {
        return this->socketPortToDatapathIdMap[socketPort];
    }
}
