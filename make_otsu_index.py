import os, sys, math, random, json

f = open("log.txt")

num_done = -1
lf = open("otsu_0.txt", "wb")
for line in f:
  fn, pset = line.strip().split(",")
  fn = fn
  sfn = "images/" + fn
  num_done += 1
  print "Working", sfn, pset
  os.system("./cotsu %s > junk_0.txt" % (sfn))
  s = open("junk_0.txt").read()
  qq = {}
  qq["file"] = fn
  qq["set"] = pset
  qq["otsu"] = s
  lf.write(json.dumps(qq) + "\n")
os.system("rm junk_0.txt")
