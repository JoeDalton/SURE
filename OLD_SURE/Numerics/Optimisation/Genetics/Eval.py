from Utils.I_O import ErrorMsg

fileName = 'Numerics/Optimisation/Genetics/Eval.py'

class MultiEval(): # Lower is better !
  xGrade  = []
  xCrit   = [] # xCrit must have one less elemnt than xGrade: n constraints and one parameter to be optimized
  xSign   = [] # 1 for minimization, -1 for maximization. Default is [1]*len(xGrade)

  def __init__(self, xGrade, xCrit, xSign=[]):
    self.xGrade     = xGrade
    self.xCrit      = xCrit
    if xSign == []:
      self.xSign = [1]*len(xGrade)
    else:
      self.xSign = xSign
    if len(xCrit)>=len(xGrade):
      ErrorMsg("xCrit too long", fileName)
    elif len(xCrit)<len(xGrade)-1:
      ErrorMsg("xCrit too short", fileName)

    
  def __eq__(self, other):
    # Fool proofing : check lengths of lists and that xSign are the same
    i = 0
    while i < len(self.xCrit): # Check for all constraints first
      if self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        i += 1
      elif self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] > other.xCrit[i]*other.xSign[i]:
        return False
      elif self.xGrade[i]*self.xSign[i] > self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        return False
      else:
        if self.xGrade[i]*self.xSign[i] == other.xGrade[i]*other.xSign[i]:
          return True
        else:
          return False

    # We get here only if all the constraints are satisfied
    if self.xGrade[-1]*self.xSign[-1] == other.xGrade[-1]*other.xSign[-1]:
      return True
    else:
      return False


  def __lt__(self, other):
    # Fool proofing : check lengths of lists and that xSign are the same
    i = 0
    while i < len(self.xCrit): # Check for all constraints first
      if self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        i += 1
      elif self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] > other.xCrit[i]*other.xSign[i]:
        return True
      elif self.xGrade[i]*self.xSign[i] > self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        return False
      else:
        if self.xGrade[i]*self.xSign[i] < other.xGrade[i]*other.xSign[i]:
          return True
        else:
          return False

    # We get here only if all the constraints are satisfied
    if self.xGrade[-1]*self.xSign[-1] < other.xGrade[-1]*other.xSign[-1]:
      return True
    else:
      return False


  def __gt__(self, other):
    # Fool proofing : check lengths of lists and that xSign are the same
    i = 0
    while i < len(self.xCrit): # Check for all constraints first
      if self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        i += 1
      elif self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] > other.xCrit[i]*other.xSign[i]:
        return False
      elif self.xGrade[i]*self.xSign[i] > self.xCrit[i]*self.xSign[i] and \
            other.xGrade[i]*other.xSign[i] <= other.xCrit[i]*other.xSign[i]:
        return True
      else:
        if self.xGrade[i]*self.xSign[i] > other.xGrade[i]*other.xSign[i]:
          return True
        else:
          return False

    # We get here only if all the constraints are satisfied
    if self.xGrade[-1]*self.xSign[-1] > other.xGrade[-1]*other.xSign[-1]:
      return True
    else:
      return False


  def __le__(self, other):
    return not other < self

  def __ge__(self, other):
    return not other > self

  def Print(self):
    myStr = ""
    for i in range(len(self.xCrit)):
      if self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i]:
        myStr += "T"
      else:
        myStr += "F"
    myStr += "_" + str(self.xGrade[-1])
    print(myStr)
    print(self.xGrade)
    print(self.xSign)
    print(self.xCrit)

  def Display(self):
    myStr = ""
    for i in range(len(self.xCrit)):
      if self.xGrade[i]*self.xSign[i] <= self.xCrit[i]*self.xSign[i]:
        myStr += "T"
      else:
        myStr += "F"
    myStr += "_" + str(self.xGrade[-1])
    return myStr
