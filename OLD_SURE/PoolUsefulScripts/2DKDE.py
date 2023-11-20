import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

# Create random data
np.random.seed(0)
data = np.random.multivariate_normal((0, 0), [[0.8, 0.05], [0.05, 0.7]], 100)
x = data[:, 0]
y = data[:, 1]
xmin, xmax = -3, 3
ymin, ymax = -3, 3

# Peform the kernel density estimation
values = np.vstack([x, y])
kde = st.gaussian_kde(values)




# Find the isoValue so that the probability of being inside this iso-contour is at the specified level (similar to XX% confidence interval)
desiredProbability = 0.95

# Evaluation of the probability inside a contour is done with Monte-Carlo integration (brute-force ; there may be subtler ways to do it but it is late and I want to sleep so judgmental comments are unwelcome)
sampleSize      = 100000
samples         = kde.resample(size=sampleSize)
pdfValues       = kde(samples)
highestPDFValue = max(pdfValues)

def EstimateProbability(isoValue):
  # Filter the samples to retain only those inside your contour
  inSamples = pdfValues > isoValue #gives array of booleans
 
  # Compute integral 
  integral = inSamples.sum() / float(sampleSize)
  return integral


def FindOptimalIsoValue(desProba, fun, highVal):
  threshold = 0.0001
  counter   = 0
  residual  = 1000000000.0
  low = 0.0
  high = highVal
  while abs(residual) > threshold:
    if counter > 1000:
      print("Not converged. Reaching this state is highly unlikely, so I recommend you play the lottery if it happens.")
      break
    counter += 1
    guess = (high + low) / 2.0
    residual = fun(guess) - desProba
    if residual > 0.0: # probability inside contour is higher than desired -> must tighten the contour -> must increase the isoValue
      low   = guess
    else:
      high  = guess 
  return guess


myIsoValue = FindOptimalIsoValue(desiredProbability, EstimateProbability, highestPDFValue)
print("Iso value  = " + str(myIsoValue))
print("Proba      = " + str(EstimateProbability(myIsoValue)))


# Plot 2D kernel density estimation with the isoline computed above
xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
positions = np.vstack([xx.ravel(), yy.ravel()])
f = np.reshape(kde(positions).T, xx.shape)
fig = plt.figure()
ax = fig.gca()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
# Contourf plot
cfset = ax.contourf(xx, yy, f, cmap='Blues')
#cset = ax.contour(xx, yy, f, colors='k')
# Or kernel density estimate plot instead of the contourf plot
#ax.imshow(np.rot90(f), cmap='Blues', extent=[xmin, xmax, ymin, ymax])
# Contour plot
cset = ax.contour(xx, yy, f, levels=[myIsoValue], colors='r')
# Label plot
def fmt(x):
  return str(desiredProbability * 100.0) + "%"

ax.clabel(cset, inline=1, fmt=fmt, fontsize=10)
ax.set_xlabel('Y1')
ax.set_ylabel('Y0')
plt.show()
