from scipy import integrate
#This integrates laplace polynomials to check orthoganalness

class Memoize: #Use this to cashe results of Lpoly to speed it up A LOT. 
    def __init__(self, f): #In fact, without this the integrals take impractially long
        self.f = f           #This is (almost) ripped straight off of stack exchange.
        self.memo = {}       #But essentially this will save results so we don't have to
    def __call__(self, *args): #Keep computing recursive levels we've already computed. 
        if not args in self.memo:
            self.memo[args] = self.f(*args)
        return self.memo[args]


def Lpol(x,n): #Blatant rip off of the terms in the notes.
	if n == 0:
		return 1
	if n == 1:
		return x
	term1 = (2*n -1) * x * Lpol(x,n-1)
	term2 = (n-1) * Lpol(x,n-2)
	return ((term1-term2)/n)

Lpol = Memoize(Lpol) #Remember things that it calculated
f=open("poly2.txt",'w') #Output data for a pretty plot
for i in range(0,200):
	f.write("{} {}\n".format(i/100-1,Lpol(i/100-1,10)))


print(Lpol(.5,4))
file = open("errors.txt",'w')
for j in range(0,100,10):
    for i in range(0,100,1):
        if i == j:
            continue
        function = lambda x: Lpol(x,j) * Lpol(x,i)
        ans = integrate.romberg(function,-1,1.0)
        file.write("{}     {}     {}\n".format(j,i,abs(ans))) #output data for plotting. 
    file.write("\n")
    file.write("\n")
	#print(j)
		#print(ans,i,i+1)

for i in range(0,99,1): #this shows zero inner product with successive polynomials
    function = lambda x: Lpol(x,i) * Lpol(x,i+1)
    ans = integrate.romberg(function,-1,1.0)
    print(ans)
