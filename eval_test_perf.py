import os, math, sys, random


classes = {}
for i in range(4):
  classes[i] = []

coeffs = [0.0 for i in range(630)]

for line in open("coeffs.txt"):
  a,v = line.strip().split()
  if int(a) == 0: continue
  coeffs[int(a)-1] = float(v)

test = []
for line in open("transcript.txt"):
  pieces = line.strip().split() 
  cl = int(pieces[0])
  fn = pieces[1]
  v = [float(x.strip().split(":")[-1]) for x in pieces[2:]]
  assert len(v) == len(coeffs)
  dumbscore = sum([cc*vv for cc,vv in zip(coeffs,v)])
  # TODO clean up this debugging hack
  #if fn.find("9327") > -1:
  #  for cc,vv,ss in zip(coeffs,v,pieces[2:]):
  #    print "%10.6f %10.6f %10s" % (cc,cc*vv,ss)
  #  print "==>", dumbscore
  classes[cl].append(v)
  if cl in [2,3]: test.append( (dumbscore, random.uniform(0,1), cl, fn) )

test.sort()
test.reverse()

npos = len(classes[2])
nneg = len(classes[3])

poskept = 0
negreject = nneg

targ = 0.85
for sco, ro, cl, fn in test:
  # print sco, fn  # debugging hack
  if poskept >= targ*npos: break
  if cl == 2: poskept += 1
  else: negreject -= 1
  lastscore = sco

print "        J score: %8.4f" % ((targ + 3*negreject/float(nneg))/4.0)
odds = math.exp(lastscore)
print "    Threshold %%: %8.4f" % (odds/(1+odds))
print "Threshold score: %8.4f" % lastscore
