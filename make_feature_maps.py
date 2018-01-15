import math, os, sys, random, json

fns = []
for i in [0]: #,1,2]
  fns.append ( "otsu_%i.txt" % i )

# ==== cofigurations ==========================================================

v1_logo = ("Mi.50_v1app/", ["103"], ["101"])
v1_text = ("Mi.50_v1app/", ["104"], ["102"])

v2_Mi50_logo = ("Mi.50/", ["103"], ["101"])
v2_Mi50_text = ("Mi.50/", ["104"], ["102"])

v2_Nii40_logo = ("Nii.40/", ["109"], ["101"])
v2_Nii40_text = ("Nii.40/", ["110"], ["102"])

# ==== choose a configuration =================================================

#config = v1_logo # 69, 70, 87
#config = v1_text # 57, 52, 73

#config = v2_Nii40_logo # 67, 73, 81, 86
#config = v2_Nii40_text # 65, 84, 92, 96 

config = v2_Mi50_logo # 83, 84, 89
#config = v2_Mi50_text # 75, 82, 82, 89 

# ============================================================================

prefix, posclass, negclass = config

counts  = {}

def ovec (o):
  o = o.strip().split("\n")
  buff = []
  buff2 = []
  for x in o:
    xs = x.split()
    yL = float(xs[1])
    yR = float(xs[2])
    buff.append(yL)
    buff2.append(yR)
  return buff + buff2

def vectorize(J):
  oyflash = ovec(J["otsu"])
  return oyflash 

# == static seed for reproducible test/train split ===========================
random.seed(200)
pos = []
neg = []
buff = []
num_filtered = 0
min_contrast = 0.25
for fn in fns:
  for line in open(fn):
    try:
      J = json.loads(line)
    except:
      continue
    f = J["file"]
    if not f.startswith(prefix): continue
    k = J["set"]
    if k not in counts: counts[k] = 0
    counts[k] += 1
    J["vec"] = vectorize(J)
    if J["vec"][-1] < min_contrast: continue
    buff.append ( J )
    if k in negclass: neg.append(J)
    elif k in posclass: pos.append(J)

def split (x):
  rx = [(random.uniform(0,1), xv) for xv in x]
  rx.sort()
  K = (len(rx)*2)/3
  return [xv for r,xv in rx[:K]], [xv for r,xv in rx[K:]]

pos_train, pos_test = split(pos)
neg_train, neg_test = split(neg)

# == feature statistics ======================================================

def mstd(vns, p_mu = 0.0, p_std = 1.0, p_n = 0):
  totv = 0 + p_mu*p_n
  totv2 = 0 + (p_std*p_std + p_mu*p_mu)*p_n
  totn = 0 + p_n
  if isinstance(vns[0], tuple):
    for v,n in vns:
      totv += v*n
      totv2 += v*v*n
      totn += n
  else:
    for v in vns:
      totv += v
      totv2 += v*v
      totn += 1.0
  mu = totv/float(totn)
  va = totv2/float(totn) - mu*mu
  return mu, va**0.5

def make_buckets (vals):
  m,s = mstd(vals)
  D = s*0.25
  buff = [ (m-D/2.0 +(i-30)*D,m+D/2.0+(i-30)*D) for i in range(61) ]
  return buff

def bucketize (buck, v):
  i = len(buck)/2
  L,R = buck[i]
  if L <= v <= R: return i
  D = -1 if v < L else 1
  while 0 < i < len(buck)-1:
    i += D
    L,R = buck[i]
    if L <= v <= R: return i
  return i

# == map all the relative densities ==========================================

scobuff = []
maps = {}
out_map = []
for K0 in range(36):
  buck0 = make_buckets( [x["vec"][K0] for x in pos_train+neg_train] )
  for K1 in range(K0+1,36):
    buck1 = make_buckets( [x["vec"][K1] for x in pos_train+neg_train] )
    NB0 = len(buck0)
    pos_counts = [[0 for i in range(NB0)] for j in range(NB0)]
    neg_counts = [[0 for i in range(NB0)] for j in range(NB0)]
 
    for co,x in [ (pos_counts, pos_train), (neg_counts, neg_train) ]:
      for v in x:
        vec = v["vec"]
        b1 = bucketize(buck0,vec[K0])
        b2 = bucketize(buck1,vec[K1])
        co[b1][b2] += 1
    f = open("map.txt","w")
    f.write("%i %i\n" % (len(pos_counts[0]), len(pos_counts)))
    for i in range(len(pos_counts)):
      for j in range(len(pos_counts[0])):
        n1,n2 = (pos_counts[i][j], neg_counts[i][j])
        if max(n1,n2) == 0: continue
        f.write("%i %i %i %i\n" % (i,j,n1,n2))
    f.close() 
    os.system("python map_feature.py map.txt")
    J = json.loads(open("mapdata.txt").read())
    scobuff.append( (J[-1], K0, K1) )
    scobuff.sort()
    maps[(K0,K1)] = (buck0, buck1, J)
    out_map.append ( [K0,K1,buck0,buck1,J]  )
    print "... %4i %4i"  % (K0,K1)
  print "=", K0

# == feature transcript ======================================================

f = open("transcript.txt","wb")
for x,s in zip( [pos_train, neg_train, pos_test, neg_test], [0,1,2,3]):
  for pc,p in enumerate(x):
    v = p["vec"]
    fn = p["file"]
    keyn = 0
    buff = []
    for K0 in range(36):
      for K1 in range(K0+1,36):
        buck0, buck1, J = maps[(K0,K1)]
        b0 = bucketize(buck0,v[K0])
        b1 = bucketize(buck1,v[K1])
        sco = "%1.4f" % J[b0][b1]
        key = "%s:%s:%s:%s:%s:" % (K0,b0,K1,b1,keyn)
        buff.append(key + str(sco))
        keyn += 1
    f.write(str(s) + " " + fn + " ")
    f.write(" ".join(buff) + "\n")
f.close()

open("master_map.json","wb").write(json.dumps(out_map))

