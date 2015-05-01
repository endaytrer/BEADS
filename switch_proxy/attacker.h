/******************************************************************************
* Author: Samuel Jero <sjero@purdue.edu>
* SDN Switch-Controller Proxy
******************************************************************************/
#ifndef _ATTACKER_H
#define _ATTACKER_H
#include "sw_proxy.h"
#include "half_conn.h"
#include <map>
#include <list>

#define ACTION_ID_ERR			(-1)
#define ACTION_ID_MIN			0
#define ACTION_ID_DROP			0
#define ACTION_ID_DELAY			1
#define ACTION_ID_DUP			2
#define ACTION_ID_LIE			3
#define ACTION_ID_PRINT			4
#define ACTION_ID_MAX			4


typedef std::map<int, int> amap_t;
typedef std::map<int, std::map<int, int> > aamap_t;
typedef std::map<int,std::map<int, std::map<int, int> > > aaamap_t;
typedef std::map<uint64_t, std::map<int,std::map<int, std::map<int, int> > > > aaaamap_t;
typedef std::map<int, std::map<uint64_t, std::map<int,std::map<int, std::map<int, int> > > > > aaaaamap_t;

class Attacker{
	private:
		Attacker();

	public:
		~Attacker();
		static Attacker& get();
		bool addCommand(Message m);
		pkt_info doAttack(pkt_info pk);

	private:
		int normalize_ofp_ver_str(char *v);
		int normalize_ofp_msg_type(int ver, char *s);
		int normalize_action_type(char *s);
		int normalize_cid(char *s);
		uint64_t normalize_dpid(char *s);

		pthread_rwlock_t lock;
		// <cid, <dpid, <of_version, <pkt_type, <action, ID> > > >
		aaaaamap_t actions_map;
		// <ID, list_of_parameters >
		std::map<int, std::list<int> > params;
};


#endif
