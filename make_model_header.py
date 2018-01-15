import json, sys

J = json.loads(open("master_map.json").read())

K = {}
for line in open("coeffs.txt"):
  id, w = line.strip().split()
  id = int(id)
  w = float(w)
  K[id] = w

print "float model[] = {"
for q in range(len(J)):
  jj = J[q]
  w = K[q+1]
  K0,K1,buck0,buck1,k0k1map = jj
  # 4*61
  thing = [ [x[0] for x in buck0], [x[1] for x in buck0] ]
  thing += [ [x[0] for x in buck1], [x[1] for x in buck1] ]
  for row in k0k1map:
    # 61*61
    for i in range(len(row)):
      row[i] *= w
    thing.append(row)
  print "  //", q
  for t in thing:
    for i in range(8):
      out = t[8*i:8*(i+1)]
      o1 = out[:6]
      o2 = out[6:]
      f1 = ",".join(["%6.4e"]*len(o1))
      f2 = ",".join(["%6.4e"]*len(o2))
      if len(o1) > 0:
        print "  " + f1 % tuple(o1), ","
      if len(o2) > 0:
        print "  " + f2 % tuple(o2), ","
print "};"

