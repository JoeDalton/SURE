import math
import numpy as np
from numpy.linalg import inv
from scipy.special import erfinv
from scipy.stats import linregress
from Numerics.Algebra import NRMSE, NMAE
#from sklearn.linear_model import LinearRegression
from SurrogateModelling.SurrogateModel import SurrogateModel as SM
from SurrogateModelling.SurrogateTools.Autocorellation_functions import *
from SurrogateModelling.SurrogateTools.optimisationTools import *
from SurrogateModelling.PCE import *
from Utils.I_O import prt

fileName = 'SurrogateModelling/Kriging.py'

#TODO: Refactor to unify trend models (constant, linear, PCE)

class PCK_model(SM):
    # Details: Class that defines a universal kriging surrogate model

    def __init__(self, af, x, y, distribution, PCEdeg=2, givenTheta=None, verbosity=0, **kwargs):
        # Details: Constructor for Polynomial Chaos-Kriging class
        #
        # inputs:
        # af         - Function handle for autocorrelation function (must be a string)
        # x          - Input sample values
        # y          - Output sample values
        # PCEdeg     - (optional) degree of PCE used as trend. Default is 2
        # optiStrat  - (optional) optimization strategy for hyperparameters
        # optiTarget - (optional) Objective function to minimize for the optimization of hyperparameters
        # givenTheta - (optional) Define Kriging model with specific hyperparameter values (must be a vector with dimension iDim)
        # nugget     - (optional) nugget value: variance of noisy samples with assumed normal distribution. Can be very small (1E-10) just to make R reversible
        #
        # outputs:
        # self       - Ordinary kriging class
        # m          - Number of samples
        # iDim       - Dimension of input
        if "ID" in kwargs:
            self.ID = kwargs["ID"]
        else:
            self.ID = ""

        # Set verbosity
        if "verbosity" in kwargs:
            self.verbosity = kwargs["verbosity"]  # Should be an integer
        else:
            self.verbosity = 0

        # Set optimization strategy
        if "optiStrat" in kwargs:
            self.optiStrat = kwargs["optiStrat"]
        else:
            self.optiStrat = 'default'

        # Set optimization target
        if "optiTarget" in kwargs:
            self.optiTarget = kwargs["optiTarget"]
        else:
            self.optiTarget = 'MLE'

        # Set nugget
        if "nugget" in kwargs:
            self.nugget = kwargs["nugget"]
        else:
            self.nugget = 0.05

        prt('', 'green', self.verbosity)
        prt('Building Universal Kriging ' + self.ID + '...', 'green', self.verbosity)

        if PCEdeg == 2:
            self.PCE = PCE(ID='2', verbosity=self.verbosity - 1, distribution=distribution, checkBasis=False)
            self.PCE.DefinePolyBasis(quasiNorm=1.0)
            self.PCE.DefineTruncatureStrat(strategy="Fixed", maxTotDegree=2)
        elif PCEdeg == 3:
            self.PCE = PCE(ID='3', verbosity=self.verbosity - 1, distribution=distribution, checkBasis=False)
            self.PCE.DefinePolyBasis(quasiNorm=0.5)  # Get only third order and second-order interactions
            self.PCE.DefineTruncatureStrat(strategy="Fixed", maxTotDegree=4)
        else:
            print("Too high polynomial degree")
            quit()
        self.PCE.DefineEvalStrat(evalStrategy="Regression")
        self.PCE.Build(x, y)

        xPCEEval = np.zeros_like(y)
        # Substract PCE to data
        for i in range(len(xPCEEval)):
            xPCEEval[i] = self.PCE.Evaluate(x[i, :])[0]

        if givenTheta is None:
            self.OK = OK_model(af, x, y - xPCEEval, ID='OK', optiStrat=self.optiStrat, optiTarget=self.optiTarget,
                               nugget=self.nugget, verbosity=self.verbosity - 1)
            print(self.OK.theta)
        else:
            self.OK = OK_model(af, x, y - xPCEEval, ID='OK', givenTheta=givenTheta, nugget=0.05, verbosity=verbosity)

    def Evaluate(self, x0):
        return self.predict(x0)

    def predict(self, x0, **kwargs):
        # Details: Prediction of the surrogate model for input x0
        #
        # inputs:
        # self              - Universal kriging class object
        # x0                - Input value/vector
        # conf (opt)        - confidence interval (in percentage) - Warning, if this option is chosen, x0 must be a single sample
        #
        # outputs:
        # mean              - Mean predition output
        # variance          - Variance prediction output
        # conf_l            - Lower prediction at given confidence
        # conf_u            - Upper prediction at given confidence

        if "conf" in kwargs:
            pceVal = self.PCE.Evaluate(x0)[0]
            OKmean, OKvar, OKl, OKu = self.OK.Evaluate(x0, conf=kwargs["conf"])
            return pceVal + OKmean, OKvar, pceVal + OKl, pceVal + OKu
        else:
            return self.PCE.Evaluate(x0)[0] + self.OK.Evaluate(x0)


class OK_model(SM):

    def __init__(self, af, x, y, **kwargs):
        # Details: Constructor for ordinary Kriging class
        #
        # inputs:
        # af         - Function handle for autocorrelation function (must be a string)
        # x          - Input sample values
        # y          - Output sample values
        # optiStrat  - (optional) optimization strategy for hyperparameters
        # optiTarget - (optional) Objective function to minimize for the optimization of hyperparameters
        # givenTheta - (optional) Define ordinary Kriging model with specific hyperparameter values (must be a vector with dimension iDim)
        # nugget     - (optional) nugget value: variance of noisy samples with assumed normal distribution. Can be very small (1E-10) just to make R reversible

        # Warning : default optimisation may be unstable when used with NRMSE. Please use "GeomDichotomy" instead with NRMSE if necessary. Only possible for 1D situations though...

        # Set ID
        if "ID" in kwargs:
            self.ID = kwargs["ID"]
        else:
            self.ID = ""

        # Set verbosity
        if "verbosity" in kwargs:
            self.verbosity = kwargs["verbosity"]  # Should be an integer
        else:
            self.verbosity = 0

        prt('', 'green', self.verbosity)
        prt('Building Kriging ' + self.ID + '...', 'green', self.verbosity)

        self.af = af
        if af == "Matern32":
            self.auto_correlation_function = Matern32_matrix
        elif af == "Matern52":
            self.auto_correlation_function = Matern52_matrix
        elif af == "Exp":
            self.auto_correlation_function = Ex_matrix
        elif af == "SqExp":
            self.auto_correlation_function = SqEx_matrix
        elif af == "CubicSpline":
            self.auto_correlation_function = CubicSpline_matrix
        elif af == "Linear":
            self.auto_correlation_function = Linear_matrix
        else:
            raise Exception('Unknown covariance kernel : "' + af + '"')

        self.X = x.copy()  # x = [x1 y1 z1; x2 y2 z2]
        self.Y = y.copy()
        if hasattr(x[0], '__len__'):
            self.iDim = len(x[0])
        else:
            self.iDim = 1

        self.m = np.shape(x)[0]
        self.F = np.ones(self.m).T

        # Set optimisation method as prescribed or set default
        if "optiStrat" in kwargs:
            self.theta_opti_technique = kwargs["optiStrat"]
        else:
            self.theta_opti_technique = "default"
        if "optiTarget" in kwargs:
            self.theta_target = kwargs[
                "optiTarget"]  # For now, only "MLE" (default), "NRMSE" (expensive) and "NMAE" (expensive and not recommended a priori)
        else:
            self.theta_target = "default"

        # Set nugget
        if "nugget" in kwargs:
            self.nugget = kwargs["nugget"]
        else:
            self.nugget = 0.0

        # Substract trend if relevant
        if "trend" in kwargs:
            self.trendType = kwargs["trend"]
            if self.trendType == "linear":
                # Compute linear regression
                self.trend = linregress(self.X, self.Y)
                tosubstract = self.trend.intercept + self.trend.slope * self.X
                self.Y -= tosubstract
            elif self.trendType == "none":
                pass  # regression towards the mean
            else:
                raise Exception('Unknown trend type : "' + self.trend + '"')
        else:
            self.trendType = "none"

        # Set correlation as prescribed or optimize it by default
        if "givenTheta" in kwargs:
            self.theta = kwargs["givenTheta"]
        else:
            self.theta = self.optimise_theta()

        # Compute auto-correlation function, modified with a nugget for noisy data
        self.R = self.compute_R(self.theta)
        self.R = self.R + np.identity(self.m) * self.nugget

        # Compute a priori mean and variance
        self.beta_hat = self.compute_beta_hat(self.R)
        self.sigma_sq_hat = self.compute_sigma_sq_hat(self.R, self.beta_hat)

        prt('Kriging ' + self.ID + ' computed', 'green', self.verbosity)

    def Evaluate(self, x0):
        return self.predict(x0)

    def predict(self, x0, **kwargs):
        # Details: Prediction of the surrogate model for input x0
        #
        # inputs:
        # self              - Ordinary kriging class object
        # x0                - Input value/vector
        # conf (opt)        - confidence interval (in percentage) - Warning, if this option is chosen, x0 must be a single sample
        #
        # outputs:
        # mu_hat            - Mean predition output
        # sigma_Y_sq_hat    - Variance prediction output
        # conf_l            - Lower prediction at given confidence
        # conf_u            - Upper prediction at given confidence

        def norminv(y):
            # Quantile function for the standard normal distribution, useful to give confidence intervals
            return 1.4142135623730951 * erfinv(2.0 * y - 1.0)

        sigma_sq_hat = self.compute_sigma_sq_hat(self.R, self.beta_hat)

        mu_hat = self.compute_mu_hat(self.R, self.beta_hat, x0, self.theta)

        if self.trendType == "linear":
            tosum = self.trend.intercept + self.trend.slope * x0
            mu_hat += tosum

        if "conf" in kwargs:
            sigma_Y_sq_hat = self.compute_sigma_Y_sq_hat(self.sigma_sq_hat, x0, self.theta, self.R)
            alpha2 = 1 - kwargs["conf"] / 100.0
            uncertainty = norminv(1.0 - (alpha2 / 2.0)) * math.sqrt(sigma_Y_sq_hat)
            conf_l = mu_hat - uncertainty
            conf_u = mu_hat + uncertainty
            return mu_hat, sigma_Y_sq_hat, conf_l, conf_u
        else:
            return mu_hat

    def Smoothen(self, factor=2):
        # Experimental: Multiplies the length scales of the kriging to reduce overfitting
        self.theta *= factor

    def compute_R(self, theta):
        # Details: Obtain autocorrelation matrix
        #
        # inputs:
        # obj - Ordinary kriging class object
        # theta - Hyperparameters
        #
        # outputs:
        # R - Autocorrelation matrix

        R = self.auto_correlation_function(self.X, self.X, theta)

        # if sum(isnan(R[:])):
        #  raise Exception('NAN values')
        return R

    def compute_beta_hat(self, R):
        # Details: Obtain a priori mean
        #
        # inputs:
        # obj - Ordinary kriging class object
        # R - Autocorrelation matrix
        #
        # outputs:
        # beta_hat - A priori mean

        shouldExit = False
        try:
            invR = inv(R)
        except:
            print("Warning, R is not inversible. Please add a small nugget to correct the issue")
            shouldExit = True
        if shouldExit:
            quit()
        temp1 = np.matmul(self.F.T, invR)
        temp2 = np.matmul(temp1, self.F)

        if not hasattr(temp2, '__len__'):
            temp3 = 1.0 / temp2
            temp4 = temp3 * self.F.T
        else:
            print(fileName + ': Internal warning: the following should be a scalar:')
            print(temp2)
            temp3 = inv(temp2)
            temp4 = np.matmul(temp3, self.F.T)
        temp5 = np.matmul(invR, self.Y)
        beta_hat = np.matmul(temp4, temp5)

        return beta_hat

    def compute_sigma_sq_hat(self, R, beta_hat):
        # Details: Obtain a priori variance
        #
        # inputs:
        # obj - Ordinary kriging class object
        # R - Autocorrelation matrix
        # beta_hat - A priori mean
        #
        # outputs:
        # sigma_sq_hat - A priori variance
        if hasattr(beta_hat, '__len__'):
            temp1 = (self.Y - np.matmul(self.F, beta_hat)).T
            temp2 = self.Y - np.matmul(self.F, beta_hat)
        else:
            temp1 = (self.Y - beta_hat * self.F).T
            temp2 = self.Y - (beta_hat * self.F)
        temp3 = np.linalg.lstsq(R, temp2, rcond=None)[0]
        temp4 = np.matmul(temp1, temp3)

        sigma_sq_hat = (1.0 / self.m) * temp4
        return sigma_sq_hat

    def compute_r0(self, theta, x0):
        # Details: Obtain r0 autocorrelation vector
        #
        # inputs:
        # self - Ordinary kriging class object
        # theta - Hyperparameters
        # x0 - Input value
        #
        # outputs:
        # r0 - r0 autocorrelation vector
        r0 = self.auto_correlation_function(self.X, x0, theta)
        return r0

    def compute_mu_hat(self, R, beta_hat, x0, theta):
        # Details: Obtain the prediction mean
        #
        # inputs:
        # self - Ordinary kriging class object
        # R - Autocorrelation matrix
        # beta_hat - A priori mean
        # x0 - Input value to obtain the mean for
        # theta - Hyperparameters
        #
        # outputs:
        # mu_hat - A posteriori mean prediction
        r0 = self.compute_r0(theta, x0)
        # temp1 = self.Y - np.matmul(self.F, beta_hat)
        if not hasattr(beta_hat, '__len__'):
            temp1 = self.Y - beta_hat * self.F
        else:
            temp1 = self.Y - np.matmul(self.F, beta_hat)
        temp2 = np.linalg.lstsq(R, temp1, rcond=None)[0]
        temp3 = np.matmul(r0.T, temp2)
        mu_hat = beta_hat + temp3
        return mu_hat

    def compute_sigma_Y_sq_hat(self, sigma_sq_hat, x0, theta, R):
        # Details: Obtain the prediction variance
        #
        # inputs:
        # self - Ordinary kriging class object
        # sigma_sq_hat - A priori variance
        # x0 - Input value to obtain variance for
        # theta - Hyperparameter value
        # R - Autocorrelation matrix
        #
        # outputs:
        # sigma_Y_sq_hat - A posteriori variance prediction

        r0 = self.compute_r0(theta, x0)
        temp1 = np.linalg.lstsq(R, r0, rcond=None)[0]
        u0 = np.matmul(self.F.T, temp1) - 1.0
        temp2 = np.linalg.lstsq(R, self.F, rcond=None)[0]
        temp3 = np.matmul(self.F.T, temp2)

        try:
            temp4_1 = inv(temp3)
        except:
            temp4_1 = 1.0 / temp3
        try:
            temp4 = np.matmul(temp4_1, u0)
        except:
            temp4 = temp4_1 * u0

        try:
            temp5 = np.matmul(u0.T, temp4)
        except:
            temp5 = u0 * temp4
        temp6 = np.matmul(r0.T, temp1)

        sigma_Y_sq_hat = sigma_sq_hat * (1.0 - temp6 + temp5)

        # To avoid errors when sigma_Y_sq_hat is negative (happens to be negative and near zero at machine precision)
        myZero = np.zeros_like(sigma_Y_sq_hat)
        sigma_Y_sq_hat = np.fmax(myZero, sigma_Y_sq_hat)

        return sigma_Y_sq_hat

    def optimise_theta(self):
        # Details: Define the optimization process for the
        # hyperparameters
        #
        # inputs:
        # self - Ordinary kriging class object
        #
        # outputs:
        # theta - Optimized hyperaparameter vector

        AA = []
        b = []
        Aeq = []
        beq = []

        n = self.iDim
        lb = np.zeros(n)
        ub = np.zeros(n)

        for k in range(n):
            ite = 1
            distance = []
            for i in range(self.m):
                for j in range(self.m):
                    if not i == j:
                        try:
                            distance.append(abs(self.X[i, k] - self.X[j, k]))
                        except IndexError:
                            distance.append(abs(self.X[i] - self.X[j]))
                        ite += 1
            max_distance = max(distance)
            min_distance = min(distance)

            lb[k] = max(min_distance, 1e-5)
            ub[k] = max_distance

        if (self.theta_target == "default") or (self.theta_target == "MLE"):
            fun = self.computeMLE
        elif self.theta_target == "NRMSE":
            fun = self.computeNRMSE
        elif self.theta_target == "NMAE":
            fun = self.computeNMAE
        else:
            print("Optimisation target not supported: " + self.theta_target)
            quit()
        try:
            theta = optimisationTools(fun, self.theta_opti_technique, AA, b, Aeq, beq, lb, ub, []).x
        except AttributeError:
            theta = optimisationTools(fun, self.theta_opti_technique, AA, b, Aeq, beq, lb, ub,
                                      [])  # Case when we had a dichotomy, perhaps not a very nice way to do it...

        if self.verbosity > 0:
            print('Objective function has the value: ' + str(fun(theta)))

        return theta

    def computeMLE(self, theta):
        # Details: Define the Maximum Likelihood estimation for the
        # hyperparameters
        #
        # inputs:
        # obj - Ordinary kriging class object
        # theta - Hyperparameter input value
        #
        # outputs:
        # Psi - Value to be optimized

        R_matrix = self.compute_R(theta)
        beta_hat = self.compute_beta_hat(R_matrix)

        sigma_sq_hat = self.compute_sigma_sq_hat(R_matrix, beta_hat)
        if np.linalg.cond(R_matrix) > 1e7:
            Psi = 1e15
        else:
            Psi = 0.5 * (self.m * math.log(sigma_sq_hat) + math.log(np.linalg.det(R_matrix)))
        return Psi

    def computeNRMSE(self, theta):
        Yhat = self.compute_LOO(theta)
        res = NRMSE(self.Y, Yhat)
        return res

    def computeNMAE(self, theta):
        Yhat = self.compute_LOO(theta)
        res = NMAE(self.Y, Yhat)
        return res

    def compute_LOO(self, theta):
        Yhat = np.zeros(self.m)
        for i in range(self.m):
            newX = np.delete(self.X, i, axis=0)
            newY = np.delete(self.Y, i, axis=0)
            subKrig = OK_model(self.af, newX, newY, givenTheta=theta, nugget=self.nugget, verbosity=-1)
            Yhat[i] = subKrig.predict(self.X[i, :])
            del subKrig, newX, newY
        return Yhat

    """
    function x_new = adaptive_sampling(obj,method,A,strategy)
        % Details: Choosing the right function to create new sample based on
        % user input
        %
        % inputs:
        % obj - Ordinary kriging class object
        % method - String defining the adaptive technique
        % A - Definition of parametric space
        % strategy - Optimization technique to be used
        %
        % outputs:
        % x_new - New found sample point
        
        addpath('adaptive_techniques')
        if strcmp(method,'SSA')
            x_new = SSA_function(obj,A,strategy);
        elseif strcmp(method,'CVVor')
            x_new = CVVor_function(obj,A);
        elseif strcmp(method,'ACE')
            x_new = ACE_function(obj,A);
        elseif strcmp(method,'MIPT')
            x_new = MIPT_function(obj,A);
        elseif strcmp(method,'LOLA')
            x_new = LOLA_function(obj,A);
        elseif strcmp(method,'AME')
            x_new = AME_function(obj,A,strategy);
        elseif strcmp(method,'MEPE')
            x_new = MEPE_function(obj,A,strategy);
        elseif strcmp(method,'MASA')
            x_new = MASA_function(obj,A);
        elseif strcmp(method,'SFVCT')
            x_new = SFVCT_function(obj,A,strategy);
        elseif strcmp(method,'WAE')
            x_new = WAE_function(obj,A,strategy);
        elseif strcmp(method,'TEAD')
            x_new =TEAD_function(obj,A);
        elseif strcmp(method,'LIP')
            x_new = LIP_function(obj,A,strategy);
        elseif strcmp(method,'EI')
            x_new = EI_function(obj,A,strategy);
        elseif strcmp(method,'EIGF')
            x_new = EIGF_function(obj,A,strategy);
            
        elseif strcmp(method,'WEI')
            x_new = WEI_function(obj,A);
        elseif strcmp(method,'MSD')
            x_new = MSD_function(obj,A);
        elseif strcmp(method,'Jin_CV')
            x_new = Jin_CV_function(obj,A,strategy);
        elseif strcmp(method,'QBC_Jackknifing')
            x_new = QBC_Jackknifing_function(obj,A);
        end
    """
