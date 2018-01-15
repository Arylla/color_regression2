import math, sys, os, random, struct, json

lines = [x for x in open(sys.argv[1])]

N1,N2 = map(int, lines[0].strip().split())
N = N1

gridp = [ [0 for i in range(N)] for j in range(N) ]
gridn = [ [0 for i in range(N)] for j in range(N) ]
ngridp = [ [0 for i in range(N)] for j in range(N) ]
ngridn = [ [0 for i in range(N)] for j in range(N) ]

buff = []

for line in lines[1:]:
  i,j, np, nn = map(int, line.strip().split())
  buff.append( (i,j, np, nn) )

totp = sum([np for i,j,np,nn in buff])
totn = sum([nn for i,j,np,nn in buff])

base_ratio = float(totp)/totn
for i,j,np,nn in buff:
  gridp[i][j] = np/float(totp)
  gridn[i][j] = nn/float(totn)
  ngridp[i][j] = np
  ngridn[i][j] = nn

def gridabsmax (grid):
  buff = -1
  for i in range(len(grid)):
    for j in range(len(grid)):
      buff = max(buff, abs(grid[i][j]))
  return buff

def smooth(grid):
  N = len(grid)
  gg = [ [0 for i in range(N)] for j in range(N) ]
  for i in range(N):
    for j in range(N):
      totw = 0
      totn = 0
      if i - 1 >= 0:
        totw += 1.0
        totn += grid[i-1][j] 
      if i + 1 < N:
        totw += 1.0
        totn += grid[i+1][j] 
      if j - 1 >= 0:
        totw += 1.0
        totn += grid[i][j-1] 
      if j + 1 < N:
        totw += 1.0
        totn += grid[i][j+1] 
      gg[i][j] = (grid[i][j] + totn/totw*0.5)/1.5
  return gg

def smoothn (grid, n=1):
  while n > 0:
    grid = smooth(grid)
    n -= 1
  return grid

gridp = smoothn(gridp,50)
gridn = smoothn(gridn,50)

def gridrat (g1, g2):
  v = [[0 for i in range(len(g1))] for j in range(len(g1))]
  for i in range(len(g1)):
    for j in range(len(g1)):
      v1 = g1[i][j]
      v2 = g2[i][j]      
      v[i][j] = math.log((v1 + 1e-6)/(v2 + 1e-6))
  return v

gr = gridrat(gridp,gridn)

def drawrat (grid, fn):
  pointmap = []
  M = gridabsmax(grid)
  totalvalue = 0
  for i in range(len(grid)):
    pointmap.append([])
    row = pointmap[-1]
    for j in range(len(grid)):
      v = grid[i][j]/max(float(M),1e-6)
      row.append(v)
      totalvalue += abs(grid[i][j])*(ngridp[i][j] + ngridn[i][j])
      fr = min(abs(v),1.0)
      if v > 0: pR,pG,pB = 1.0, 1.0-fr, 1.0-fr
      else: pR,pG,pB = 1.0-fr,1.0-fr,1.0
      pR,pG,pB = [x**(1.0/2.2) for x in [pR,pG,pB]]
      pR,pG,pB = [int(round(255*x)) for x in [pR,pG,pB]]
  #pointmap.append( totalvalue/(totp+totn+1e-6)  )
  open("mapdata.txt","wb").write(json.dumps(pointmap))

drawrat(gr,"rat.png")

