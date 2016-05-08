from random import randint

def getX(n):
	s = ""
	alph = ['a','b','c','d','e','f','g','h','i','j','k','l'];
	i = 0
	while i < n:
		s += alph[randint(0,11)]
		i+=1
	return s

n = 20
x = getX(n)
y = getX(n)

with open('x.txt','w') as f:
	f.write(x)

with open('y.txt','w') as f:
	f.write(y)	
# print x,y

n = n + 1

m = [[0 for j in xrange(n)] for i in xrange(n)]
for i in xrange(n):
	for j in xrange(1):
		m[i][j] = pow(i,2)
		m[j][i] = pow(i,3)

for i in xrange(1,n):
	for j in xrange(1,n):
		m[i][j] = min(m[i][0]+pow(j,3), m[0][j] + pow(i,2))

for i in xrange(1,n):
	m[i][1] = min(m[i][1],m[i-1][0]+ (0 if x[i-1] == y[1] else 0))

for j in xrange(1,n):
	m[1][j] = min(m[1][j],m[0][j-1]+ (0 if x[1] == y[j-1] else 0))

with open('edinput.txt','w') as f:
	for i in xrange(n):
		for j in xrange(n):
			f.write(str(i)+","+str(j)+"\t"+str(m[i][j])+"\n")
			# if i==0:
			# 	f.write(str(i)+","+str(j)+"\t"+str(pow(j,3))+"\n")
			# elif j==0:
			# 	f.write(str(i)+","+str(j)+"\t"+str(pow(i,2))+"\n")