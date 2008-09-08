#!/usr/bin/env python
# this script will parse the cubep3m output file (in it's current 2006 11 29 form)
# first do a: 
# grep Redshift -A4 cubep3m_64_768_1000Mpc_061005 | grep -v Scale | grep -v Expan > ~/ts_breakdown
# then run this bad boy
import string
bd_raw=[]
file=open('ts_breakdown','r')
bd_raw.extend( file.readlines())
i=0
z = []
for l in bd_raw:
  if l[:9] == ' Redshift':
    i += 1
    l = string.split(l)
#    print l[0]
    z.append(l[3])
#print 'there were '+str(i)+' timesteps'
t1 = []
i=0
for l in bd_raw:
  if l[:2] == ' T':
    i+=1
    l = string.split(l) 
    t1.append(l[3:7])
#    print l[3:7]
#    print 't1',t1[i-1]
t2 = []
i=0
for l in bd_raw:
  if l[:2] == '  ':
    i+=1
    l = string.split(l)
    t2.append(l[0])
#    print l[0]
#    print 't2',t2[i-1] 
#print 'there were '+str(i)+' timesteps'
for j in range(i-1):
  print z[j],t1[j][0],t1[j][1],t1[j][2],t1[j][3],t2[j] 
