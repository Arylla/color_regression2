import os, math, sys, random


classes = {}
# class 0 = pos train, 1 = neg train, 2 = pos test, 3 = neg test
for i in range(4):
  classes[i] = []

def transcript_iter ():
  for line in open("transcript.txt"):
    pieces = line.strip().split() 
    cl = int(pieces[0])
    v = [float(x.strip().split(":")[-1]) for x in pieces[2:]]
    yield cl,v

test = []
for cl,v in transcript_iter():
  NV = len(v)
  dumbscore = sum(v)
  classes[cl].append(v)

npos = len(classes[0])
nneg = len(classes[1])
base_logodds = math.log(npos/float(nneg))

base_coeff = 1.0/NV
coeffs = [base_logodds] + [base_coeff for i in range(NV)]

pavLL = -99999999999999
numdone = 0
ssize = 0.1
prev_good_coeffs = [1.0] + [base_coeff for i in range(NV)]
prev_grad = [0.0 for x in prev_good_coeffs]
print "%8s %8s %8s %7s %8s" % ("LogLik", "Prob", "Iter", "Conv.", "Speed")
while True and numdone < 1000000:
  numdone += 1 
  # calculate the log likelihood and build up the gradient
  avLL = 0.0
  grad = [0.0 for i in range(NV+1)]
  for cl,v in transcript_iter():
    if cl not in [0,1]: continue
    v = [1.0] + v 
    tot = 0.0
    curr_grad = [0.0 for x in grad]
    for i in range(len(v)):
      x = v[i]
      tot += coeffs[i]*x
      curr_grad[i] = x
    odds = math.exp(tot)
    prob = odds/(1.0 + odds)
    prob_targ = 1.0 if cl == 0 else 0.0
    prob_err = prob_targ - prob
    LL = math.log(prob) if cl == 0 else math.log(1.0-prob)
    avLL += LL/float(npos+nneg)
    for j in range(len(coeffs)):
      # I've hacked the negative sign; where did I get it wrong above?
      grad[j] -= prob_err*curr_grad[j] / float(npos+nneg)
  delta = avLL - pavLL
  print "%8.4f %8.4f %8i %8.4f %8.4f" % (avLL, math.exp(avLL), numdone, (1e-6)/delta, ssize)
  if avLL < pavLL + 1e-6:
    if ssize > 1e-6 and delta < 0:
      ssize *= 0.95
      coeffs = [x for x in prev_good_coeffs]
      for j in range(len(coeffs)):
        c = coeffs[j] - ssize*prev_grad[j]
        if c > 1.0: c = 1.0
        if c < 0.0: c = 0.0
        if j == 0: c = base_logodds
        coeffs[j] = c
    else:
      break
  else:
    pavLL = avLL
    prev_good_coeffs = [x for x in coeffs]
    for j in range(len(coeffs)):
      c = coeffs[j] - ssize*grad[j]
      if c > 1.0: c = 1.0
      if c < 0.0: c = 0.0
      if j == 0: c = base_logodds
      coeffs[j] = c
    ssize *= 1.05
    prev_grad = [x for x in grad]

f = open("coeffs.txt","wb")
for j in range(len(coeffs)):
  f.write("%s %s\n" % (j,coeffs[j]))
f.close()

