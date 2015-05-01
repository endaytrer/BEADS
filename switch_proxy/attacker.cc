/******************************************************************************
* Author: Samuel Jero <sjero@purdue.edu> and Xiangyu Bu <xb@purdue.edu>
* SDN Switch-Controller Proxy
******************************************************************************/
#include "attacker.h"
#include "csv.h"
#include "args.h"
extern "C" {
#include <loci/loci.h>
#include <loci/of_object.h>
#include <loci/loci_obj_dump.h>
}
#include <stdarg.h>
#include <map>
#include <list>
using namespace std;

#define ATTACKER_ROW_NUM_FIELDS		8
#define ATTACKER_ARGS_DELIM		'&'

#define CID_ERROR				(-2)
#define CID_ALL					(-1)
#define DPID_ERROR				(-2)
#define DPID_ALL				(-1)

#define OF_VERSION_ALL			6

#define ACTION_ALIAS_DROP		"DROP"
#define ACTION_ALIAS_DELAY		"DELAY"
#define ACTION_ALIAS_DUP		"DUP"
#define ACTION_ALIAS_LIE		"LIE"
#define ACTION_ALIAS_PRINT		"PRINT"

#define OFP_MSG_TYPE_ERR		(-2)
#define OFP_MSG_TYPE_ALL		(-1)
#define OFP_MSG_TYPE_MAX		34

pthread_mutex_t ofo_print_serialization_mutex = PTHREAD_MUTEX_INITIALIZER;

Attacker::Attacker()
{
	pthread_rwlock_init(&lock, NULL);
}
Attacker::~Attacker()
{
	pthread_rwlock_destroy(&lock);
}

Attacker& Attacker::get()
{
	static Attacker me;
	return me;
}

bool Attacker::addCommand(Message m)
{
	bool ret = true;
	size_t num_fields;
	int ofp_ver;
	int msg_type;
	int action_type;
	int cid;
	int dpid;
	char **fields;
	//arg_node_t *args, *targ;
	aaaaamap_t::iterator it1;
	aaaamap_t::iterator it2;
	aaamap_t::iterator it3;
	aamap_t::iterator it4;
	amap_t::iterator it5;

	/* Parse CSV */
	fields = csv_parse(m.buff, m.len, &num_fields);
	/*for (i = 0; i < num_fields; ++i) {
		//csv_unescape(fields[i]);
		fprintf(stderr, "%d: \"%s\"\n", i, fields[i]);
	}*/

	if (num_fields != ATTACKER_ROW_NUM_FIELDS) {
		dbgprintf(0,"Adding Command: csv field count mismatch (%lu / %d).\n", num_fields, ATTACKER_ROW_NUM_FIELDS);
		ret = false;
		goto out;
	}

	if ((cid = normalize_cid(fields[1])) == CID_ERROR) {
		dbgprintf(0,"Adding Command: unrecognized CID\"%s\".\n", fields[1]);
		ret = false;
		goto out;
	}

	if ((dpid = normalize_cid(fields[2])) == DPID_ERROR) {
		dbgprintf(0,"Adding Command: unrecognized DPID \"%s\".\n", fields[2]);
		ret = false;
		goto out;
	}

	if ((ofp_ver = normalize_ofp_ver_str(fields[3])) == OF_VERSION_UNKNOWN) {
		dbgprintf(0,"Adding Command: unrecognized OFP version field \"%s\".\n", fields[3]);
		ret = false;
		goto out;
	}
	if ((msg_type = normalize_ofp_msg_type(ofp_ver, fields[4])) == OFP_MSG_TYPE_ERR) {
		dbgprintf(0,"Adding Command: unrecognized message type \"%s\" for OFP version %d.\n", fields[4], ofp_ver);
		ret = false;
		goto out;
	}

	// TODO: parse field column

	if ((action_type = normalize_action_type(fields[6])) == ACTION_ID_ERR) {
		dbgprintf(0,"Adding Command: unsupported malicious action \"%s\".\n", fields[6]);
		ret = false;
		goto out;
	}

	/* Add command */
	pthread_rwlock_wrlock(&lock);

	/* CID */
	it1 = actions_map.find(cid);
	if (it1 == actions_map.end()) {
		actions_map[cid] = aaaamap_t();
		it1 = actions_map.find(cid);
	}

	/* DPID */
	it2 = it1->second.find(dpid);
	if (it2 == it1->second.end()) {
		it1->second[dpid] = aaamap_t();
		it2 = it1->second.find(dpid);
	}

	/* OF version */
	it3 = it2->second.find(ofp_ver);
	if (it3 == it2->second.end()) {
		it2->second[ofp_ver] = aamap_t();
		it3 = it2->second.find(ofp_ver);
	}

	/* Packet Type */
	it4 = it3->second.find(msg_type);
	if (it4 == it3->second.end()) {
		it3->second[msg_type] = amap_t();
		it4 = it3->second.find(msg_type);
	}

	/* Action */
	it5 = it4->second.find(action_type);
	if (it5 == it4->second.end()) {
		it4->second[action_type] = 1;
		it5 = it4->second.find(action_type);
	} else {
		it5->second = 1;
	}

	pthread_rwlock_unlock(&lock);

out:
	csv_free(fields);
	return ret;
}

int Attacker::normalize_ofp_ver_str(char *v)
{
	int ret;

	if (strchr(v, '.')) {
		// string of format x.y or x.y.z
		if (!strncmp(v, "1.0", 3)) return OF_VERSION_1_0;
		if (!strncmp(v, "1.1", 3)) return OF_VERSION_1_1;
		if (!strncmp(v, "1.2", 3)) return OF_VERSION_1_2;
		if (!strncmp(v, "1.3", 3)) return OF_VERSION_1_3;
		if (!strncmp(v, "1.4", 3)) return OF_VERSION_1_4;
	}

	// wildcard expression
	if (v[0] == '*') return OF_VERSION_ALL;

	// the last possibility is that v is the version value, not string.
	ret = atoi(v);
	if (ret < OF_VERSION_1_0 || ret > OF_VERSION_ALL) return OF_VERSION_UNKNOWN;
	return ret;
}

int Attacker::normalize_ofp_msg_type(int ver, char *s)
{
	int ret, i;
	const char *alias;

	if (s[0] == '*') return OFP_MSG_TYPE_ALL;

	if (is_int(s)) {
		ret = atoi(s);
		if (ret > OF_OBJECT_COUNT || ret < 0)
			return OFP_MSG_TYPE_ERR;
		return ret;
	}

	for (i = 0; i < OF_OBJECT_COUNT; ++i) {
		alias = of_object_id_str[i];
		if (alias && !strcmp(alias, s)) {
			return i;
		}
	}

	// alias not found
	return OFP_MSG_TYPE_ERR;
}

int Attacker::normalize_action_type(char *s)
{
	int ret;
	if (is_int(s)) {
		ret = atoi(s);
		if (ret < ACTION_ID_MIN || ret > ACTION_ID_MAX) return ACTION_ID_ERR;
		return ret;
	}
	if (!strcmp(ACTION_ALIAS_DROP, s)) return ACTION_ID_DROP;
	if (!strcmp(ACTION_ALIAS_DELAY, s)) return ACTION_ID_DELAY;
	if (!strcmp(ACTION_ALIAS_DUP, s)) return ACTION_ID_DUP;
	if (!strcmp(ACTION_ALIAS_LIE, s)) return ACTION_ID_LIE;
	if (!strcmp(ACTION_ALIAS_PRINT, s)) return ACTION_ID_PRINT;
	return ACTION_ID_ERR;
}

int Attacker::normalize_cid(char *s)
{
	int cid = 0;

	if (s[0] == '*') return CID_ALL;

	cid = atoi(s);
	if (cid == 0) {
		return CID_ERROR;
	}
	return cid;
}

uint64_t Attacker::normalize_dpid(char *s)
{
	uint64_t dpid = 0;

	if (s[0] == '*') return DPID_ALL;

	dpid = strtoll(s,NULL,10);
	if (dpid == 0) {
		return DPID_ERROR;
	}
	return dpid;
}

int writer(void *cookie, const char *fmt, ...)
{
    va_list args;
	va_start(args, fmt);
	int ret = vfprintf(stderr, fmt, args);
	va_end(args);
	return ret;
}



pkt_info Attacker::doAttack(pkt_info pk)
{
	of_object_t* ofo = pk.ofo;
	aaaaamap_t::iterator it1;
	aaaamap_t::iterator it2;
	aaamap_t::iterator it3;
	aamap_t::iterator it4;
	amap_t::iterator it5;
	//int param;


	pthread_rwlock_rdlock(&lock);

	/* Check cid */
	it1 = actions_map.find(CID_ALL);
	if (it1 == actions_map.end()) {
		it1 = actions_map.find(pk.cid);
		if (it1 == actions_map.end()) {
			goto out;
		}
	}

	/* Check dpid */
	it2 = it1->second.find(DPID_ALL);
	if (it2 == it1->second.end()) {
		it2 = it1->second.find(pk.dpid);
		if (it2 == it1->second.end()) {
			goto out;
		}
	}

	/* Check OpenFlow Version */
	it3 = it2->second.find(pk.ofo->version);
	if (it3 == it2->second.end()) {
		goto out;
	}

	/* Check Packet Type */
	it4 = it3->second.find(pk.ofo->object_id);
	if (it4 == it3->second.end()) {
		goto out;
	}

	/* Drop */
	it5 = it4->second.find(ACTION_ID_DROP);
	if (it5 != it4->second.end()) {
		//param = it4->second[ACTION_ID_DROP];
		of_object_delete(pk.ofo);
		pk.ofo = NULL;
		goto out;
	}

	/* Debug Printing */
	it5 = it4->second.find(ACTION_ID_PRINT);
	if (it5 != it4->second.end() || sw_proxy_debug > 2) {
		pthread_mutex_lock(&ofo_print_serialization_mutex);
		if (pk.dir == STOC) {
			dbgprintf(0,"##################\nGot Message (s %llu -> c %i)\n", pk.dpid, pk.cid);
		}else {
			dbgprintf(0,"##################\nGot Message (c %i -> s %llu)\n", pk.cid, pk.dpid);
		}
		of_object_dump((loci_writer_f)&writer,NULL,ofo);
		dbgprintf(0,"##################\n");
		pthread_mutex_unlock(&ofo_print_serialization_mutex);
	}

out:
	pthread_rwlock_unlock(&lock);
	return pk;
}
