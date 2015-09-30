#!/bin/env python
# Samuel Jero <sjero@purdue.edu>
# Top-level testing script
import os
import sys
import time
from datetime import datetime
import argparse
import socket

system_home = os.path.dirname(os.path.realpath(__file__))
scripts_path = os.path.abspath(os.path.join(system_home, 'scripts'))
config_path = os.path.abspath(os.path.join(system_home, 'config'))
sys.path.insert(1,scripts_path)
sys.path.insert(0,config_path)
from scripts.test import *
import config


def main(args):
	vms_per_instance = config.controllers_per_instance + 1
	standalone = True

	#Parse args
	argp = argparse.ArgumentParser(description='Test Executor')
	argp.add_argument('-p','--port', type=int, default=config.coordinator_port)
	argp.add_argument('-c','--coordinator', type=str)
	argp.add_argument('-i','--instance', type=int, default=0)
	args = vars(argp.parse_args(args[1:]))
	instance = args['instance']
	if args['coordinator'] is not None:
		standalone = False

	print "Running Instance " + str(instance) + "..."

	#Open Log file
	lg = open(config.logs_loc.format(instance=instance), "w")
	lg.write(str(datetime.today()) + "\n")
	lg.write("Instance: " + str(instance) + "\n")

	#Determine VMs
	mininet = [instance*vms_per_instance + 1]
	controllers = list()
	for i in range(1,vms_per_instance):
		controllers.append(mininet[0] + i)
	lg.write("Mininet: " + str(mininet) + "\n")
	lg.write("Controllers: " + str(controllers) + "\n")
	lg.flush()

	#Start VMs
	print "Starting VMs..."
	#startVms(mininet, controllers)

	#Do Tests
	if standalone:
		standalone_tests(mininet,controllers,lg)
	else:
		coordinated_tests(mininet,controllers,instance, lg, (args['coordinator'], args['port']))

	#Stop VMs
	print "Stopping VMs..."
	#stopVms(mininet, controllers)

	#Close log
	lg.write(str(datetime.today()) + "\n")
	lg.close()

def reconnect(addr):
	sock = null
	while True:
		try:
			sock = socket.create_connection(addr)
			break
		except Exception as e:
			print "Failed to connect to coordinator (%s:%d): %s...retrying" % (addr[0], addr[1], e)
			time.sleep(1)
			continue
	return sock

#Comments on controller<->executor message format:
# Messages are arbitrary strings ending in \n
# use sock.send() to send and sock.makefile.readline() to receive a full message.
# readline() properly handles waiting for a full message before delivering it. On
# error, an empty string is returned.
# Arbitrary python can be passed back and forth with str=repr(data) and data=eval(str)
def coordinated_tests(mininet, controllers, instance, lg, addr):
	num = 1

	#Connect
	try:
		sock = socket.create_connection(addr)
	except Exception as e:
		print "Failed to connect to coordinator (%s:%d): %s" % (addr[0], addr[1], e)
		return
	rf = sock.makefile()

	#Loop Testing Strategies
	while True:
		#Ask for Next Strategy
		try:
			sock.send("READY %s:%d\n" %(socket.gethostname(),instance))
		except Exception as e:
			print "Failed to send on socket..."

		#Get Strategy
		line = rf.readline()
		if line=="":
			rf.close()
			sock.close()
			sock = reconnect(addr)
			rf = sock.makefile()
		strat = eval(line)

		#Test
		print strat
		doTest(mininet,controllers,strat[0],strat[1], num, lg)
		num+=1

	#Cleanup
	rf.close()
	sock.close()

def standalone_tests(mininet, controllers, lg):
	print "Starting Tests..."
	print "Test 1   " + str(datetime.today())
	res = doTest(mininet, controllers, "/root/test2.py {controllers}", ["*,*,*,*,*,CLEAR,*"], 1, lg)
	print "Test Result: " + str(res)
	print "******"
	print "Test 2   " + str(datetime.today())
	res = doTest(mininet, controllers, "/root/test2.py {controllers}", ["{controllers[0]},3,*,of_packet_in,12,CLIE,mfield=12&mval=2&act==&val=1"], 2, lg)
	print "Test Result: " + str(res)
	print "******"
	print "Test 3   " + str(datetime.today())
	res = doTest(mininet, controllers, "/root/test1.py {controllers}", ["*,*,*,*,*,CLEAR,*"], 3, lg)
	print "Test Result: " + str(res)
	print "******"
	print "Test 4   " + str(datetime.today())
	res = doTest(mininet, controllers, "/root/test1.py {controllers}", ["{controllers[0]},3,*,of_packet_in,12,CDIVERT,mfield=12&mval=3&p=100&sw=2&ctl={controllers[0]}","{controllers[0]},2,*,of_packet_in,12,CDIVERT,mfield=12&mval=3&p=100&sw=3&ctl={controllers[0]}"], 4, lg)
	print "Test Result: " + str(res)
	print "******"
	print "Test 5   " + str(datetime.today())
	res = doTest(mininet, controllers, "/root/test1.py {controllers}", ["{controllers[0]},3,*,of_packet_out,7.1.1,CDIVERT,mfield=7.1.1&mval=2&p=100&sw=1&ctl={controllers[0]}"], 5, lg)
	print "Test Result: " + str(res)
	print "******"


if __name__ == "__main__":
	main(sys.argv)
