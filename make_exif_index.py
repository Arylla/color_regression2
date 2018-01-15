import os, sys, time, math, random, json

#images = [x.strip() for x in open("images/image_list.txt")]
images = [x.strip() for x in open(sys.argv[1])]

lf = open("log.txt","wb")
num_done = -1
t0 = -1
for img_fn in images:
  num_done += 1 
  if img_fn.find("(") >= 0: continue
  if img_fn.find(" ") >= 0: continue
  if img_fn.find("egacy") > -1: continue

  bad = os.system("exiftool images/%s > junk.txt" % img_fn)
  M = {}
  M["file"] = img_fn
  for line in open("junk.txt"):
    pos = line.find(":")
    L,R = line[:pos], line[pos+1:]
    M[L.strip()] = R.strip()
  if "User Comment" not in M: continue

  # parse the user comment 
  uc = M["User Comment"]
  pieces = uc.split(",")
  d  = {}
  for p in pieces:
    pp = p.split(":")
    if len(pp) != 2: continue
    L,R = p.split(":")
    d[L] = R.strip()

  # must be flash
  if d["Count"] != "30": continue

  q = {}
  pb = d["PhotoBunch"]
  ps = d["PhotoSet"]
  q["set"] = ps
  q["file"] = M["file"]
 
  lf.write("%s,%s\n" % (q["file"],ps))
  lf.flush()
  if time.time() > t0 + 10.0:
    print "%8.2f%% " % (100*num_done/float(len(images)))
    t0 = time.time()
os.system("rm junk.txt")

