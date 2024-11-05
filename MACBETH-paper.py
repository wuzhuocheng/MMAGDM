
import math

class MACBETH:
    def __init__(self, NormalizeMatrix, Weight, q):
        '''
        :param NormalizeMatrix: 矩阵
        :param Weight: 权重
        :param Gangliang: 属性的成本或效益型
        '''
        self.Matrix = NormalizeMatrix
        self.Weight = Weight
        self.lines = len(NormalizeMatrix)  # 行数
        self.colums = len(NormalizeMatrix[0])  # 列数
        self.Gangliang = [1 for i in range(self.colums)]
        self.q = q

    def getColumsMin(self):
        Cmin = [[[0, 0], [0, 0]] for i in range(0, self.colums)]
        temp = [float('inf') for i in range(0, self.colums)]
        for i in range(0, self.lines):
            for j in range(0, self.colums):
                if temp[j] > ScoreFunction(self.Matrix[i][j], self.q):
                    temp[j] = ScoreFunction(self.Matrix[i][j], self.q)
                    Cmin[j] = self.Matrix[i][j]
        return Cmin

    def getColumsMax(self):
        Cmax = [[[0, 0], [0, 0]] for i in range(0, self.colums)]
        temp = [0 for i in range(0, self.colums)]
        for i in range(0, self.lines):
            for j in range(0, self.colums):
                if temp[j] < ScoreFunction(self.Matrix[i][j], self.q):
                    temp[j] = ScoreFunction(self.Matrix[i][j], self.q)
                    Cmax[j] = self.Matrix[i][j]
        return Cmax

    def getReferenceLevels(self):
        Cmin = self.getColumsMin()
        Cmax = self.getColumsMax()
        RL = [Cmin, Cmax]
        return RL

    def getMACBETHScores(self):
        RL = self.getReferenceLevels()
        V = [[[[0, 0], [0, 0]] for i in range(0, self.colums)] for j in range(0, self.lines)]
        for i in range(0, self.lines):
            for j in range(0, self.colums):
                temp1 = Subtract(self.Matrix[i][j], RL[0][j], self.q)
                temp2 = Subtract(RL[1][j], RL[0][j], self.q)
                temp3 = Divide(temp1, temp2, self.q)
                V[i][j] = CountMultiplication(temp3, 5, self.q)
        return V

    def getMACBETHWeightScores(self):
        V = self.getMACBETHScores()
        VW = [[[[0, 0], [0, 0]] for i in range(0, self.colums)] for j in range(0, self.lines)]
        for i in range(0, self.lines):
            for j in range(0, self.colums):
                VW[i][j] = CountMultiplication(V[i][j], self.Weight[j], self.q)
        return VW

    def getOverallScore(self):
        VW = self.getMACBETHWeightScores()
        OS = [[[0, 0], [0, 0]] for i in range(0, self.lines)]
        for i in range(0, self.lines):
            temp = VW[i][0]
            for j in range(1, self.colums):
                temp = Add(temp, VW[i][j], self.q)
            OS[i] = temp
        return OS

    def getResult(self):
        OS = self.getOverallScore()
        SF = [0 for i in range(0, self.lines)]
        for i in range(0, self.lines):
            SF[i] = round(ScoreFunction(OS[i], self.q), 6)
        Ranking = []
        temp = SF.copy()
        temp.sort(reverse=True)
        for i in range(0, self.lines):
            for j in range(0, self.lines):
                if temp[i] == SF[j]:
                    Ranking.append('A' + str(j + 1))
        SF=[i/sum(SF) for i in SF]
        return SF




def Add(Fnum1, Fnum2, q):
    """
    区间值广义正交模糊数加法运算
    :param Fnum1:
    :param Fnum2:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = pow(pow(Fnum1[0][0], q) + pow(Fnum2[0][0], q) - (pow(Fnum1[0][0], q) * pow(Fnum2[0][0], q)), 1 / q)
    temp[0][1] = pow(pow(Fnum1[0][1], q) + pow(Fnum2[0][1], q) - (pow(Fnum1[0][1], q) * pow(Fnum2[0][1], q)), 1 / q)
    temp[1][0] = Fnum1[1][0] * Fnum2[1][0]
    temp[1][1] = Fnum1[1][1] * Fnum2[1][1]
    return temp

def Subtract(Fnum1, Fnum2, q):
    """
    区间值广义正交模糊数减法运算
    :param Fnum1:
    :param Fnum2:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = Fnum1[0][0] * Fnum2[1][0]
    temp[0][1] = Fnum1[0][1] * Fnum2[1][1]
    temp[1][0] = pow(pow(Fnum1[1][0], q) + pow(Fnum2[0][0], q) - pow(Fnum1[1][0], q) * pow(Fnum2[0][0], q), 1 / q)
    temp[1][1] = pow(pow(Fnum1[1][1], q) + pow(Fnum2[0][1], q) - pow(Fnum1[1][1], q) * pow(Fnum2[0][1], q), 1 / q)
    return temp

def Multiply(Fnum1, Fnum2, q):
    """
    区间值广义正交模糊数乘法运算
    :param Fnum1:
    :param Fnum2:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = Fnum1[0][0] * Fnum2[0][0]
    temp[0][1] = Fnum1[0][1] * Fnum2[0][1]
    temp[1][0] = pow(pow(Fnum1[1][0], q) + pow(Fnum2[1][0], q) - pow(Fnum1[1][0], q) * pow(Fnum2[1][0], q), 1 / q)
    temp[1][1] = pow(pow(Fnum1[1][1], q) + pow(Fnum2[1][1], q) - pow(Fnum1[1][1], q) * pow(Fnum2[1][1], q), 1 / q)
    return temp


def Divide(Fnum1, Fnum2, q):
    """
    区间值广义正交模糊数除法运算
    :param Fnum1:
    :param Fnum2:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = pow(pow(Fnum1[0][0], q) + pow(Fnum2[1][0], q) - (pow(Fnum1[0][0], q) * pow(Fnum2[1][0], q)), 1 / q)
    temp[0][1] = pow(pow(Fnum1[0][1], q) + pow(Fnum2[1][1], q) - (pow(Fnum1[0][1], q) * pow(Fnum2[1][1], q)), 1 / q)
    temp[1][0] = Fnum1[1][0] * Fnum2[0][0]
    temp[1][1] = Fnum1[1][1] * Fnum2[0][1]
    return temp


def CountMultiplication(Fnum, a, q):
    """
    区间值广义正交模糊数数乘运算
    :param Fnum:
    :param a:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = pow(1 - pow(1 - pow(Fnum[0][0], q), a), 1 / q)
    temp[0][1] = pow(1 - pow(1 - pow(Fnum[0][1], q), a), 1 / q)
    temp[1][0] = pow(Fnum[1][0], a)
    temp[1][1] = pow(Fnum[1][1], a)
    return temp


def Power(Fnum, a, q):
    """
    区间值广义正交模糊数幂运算
    :param Fnum:
    :param a:
    :param q:
    :return:
    """
    temp = [[0, 0], [0, 0]]
    temp[0][0] = pow(Fnum[0][0], a)
    temp[0][1] = pow(Fnum[0][1], a)
    temp[1][0] = pow(1 - pow((1 - pow(Fnum[1][0], q)), a), q)
    temp[1][1] = pow(1 - pow((1 - pow(Fnum[1][1], q)), a), q)
    return temp


def ScoreFunction(Fnum, q):
    """
    得分函数
    :param Fnum:
    :param q:
    :return:
    """
    a = Fnum[0][0]**q
    b = Fnum[0][1]**q
    c = Fnum[1][0]**q
    d = Fnum[1][1]**q
    temp = (((a+b+c+d)/2)+math.exp(((a-c)+(b-d)))-math.exp((b-a)-(d-c))) / math.exp(2)
    # temp = (a+b-c-d) / 2
    # temp = 1/4*((1/2*(a+b+c+d))**2+(a+b-c-d)+(1-(b-a))**2+(d-c)**2)
    return temp


def ExactValueFunction(Fnum, q):
    """
    精确值函数
    :param Fnum:
    :param q:
    :return:
    """
    temp = pow(q + 10, (1 / 2) * (Fnum[0][0] + Fnum[0][1] + Fnum[1][0] + Fnum[1][1]))
    return temp



