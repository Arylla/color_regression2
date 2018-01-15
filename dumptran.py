s = "9327"
for line in open("transcript.txt"):
  if line.find(s) < 0: continue
  print line.strip().split()[1]
  x = line.strip().split()[2:]
  for y in x:
    print y
  break
