#!/usr/bin/python
import re
import sys
maxiter = 30
samples = 100
rankSize = 10
local_n = samples / rankSize
pat1 = re.compile(r"data rss: (.*)\n")
pat2 = re.compile(r"velset.*l1norm: (.*),.*l2norm: (.*)\n")
pat3 = re.compile(r"print HA\' \* U \* SSqInv \* U' \* \(D - HA\).*?:.*?:.*?:.*?:.*?:.*?:(.*?)[0-9]+-[0-9]+.*?DEBUG", re.S)
res1 = []
res2 = []
res3 = []
for rank in range(0,rankSize):
	fname = "enfwi-damp-" + str(rank) + ".log"
	f = open(fname)
	content =  f.read()
	res1.append(pat1.findall(content))
	res2.append(pat2.findall(content))
	res3.append(pat3.findall(content))
	f.close()

for i in range(0,maxiter):
	for rank in range(0,rankSize):
		for j in range(0,local_n):
			pos = i * local_n + j
			print res2[rank][pos][0],' ',
	print
	for rank in range(0,rankSize):
		for j in range(0,local_n):
			pos = i * local_n + j
			print res2[rank][pos][1],' ',
	print
	for rank in range(0,rankSize):
		for j in range(0,local_n):
			pos = i * local_n + j
			print res1[rank][pos],' ',
	print
	for rank in range(0,rankSize):
		print res3[rank][i + 1]
