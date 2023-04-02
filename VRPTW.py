import pyomo.environ as pyo
import numpy as np
import scipy.spatial.distance as sp
import matplotlib.pyplot as plt
class Instance:
    def __init__(self,vehicles,capacity,index,demand,coords,time_windows,service_time,dismatrix):
        self.vehicles = vehicles
        self.capacity = capacity
        self.index = index
        self.demand = demand
        self.coords = coords
        self.timeWindows = time_windows
        self.serviceTime = service_time
        self.dismatrix = dismatrix

    @staticmethod
    def create_instance(file_path):
        with open(file_path) as input_file:
            lines = input_file.readlines()
            vehicles, capacity = [int(i) for i in lines[4].split()]
            data = [[int(i) for i in line.split()] for line in lines[9:]]
        index = []
        demand = []
        coords = []
        timewindows = []
        servicetime = []
        for line in data:
            if line:
                index.append(line[0])
                coords.append(np.array(line[1:3]))
                demand.append(line[3])
                timewindows.append(line[4:6])
                servicetime.append(line[6])
        dismatrix = sp.squareform(sp.pdist(coords))
        return Instance(vehicles, capacity,index,demand,coords,timewindows,servicetime,dismatrix)

class VRPTW:
    def __init__(self,instance):
        self.instance = instance
        self.buildmodel()
        #arcs = self.solve()
        #self.draw(arcs)
    def buildmodel(self):
        #create model
        self.model = pyo.ConcreteModel("VRPTW")
        #create set
        self.model.N = pyo.Set(initialize=self.instance.index)
        self.model.C = self.model.N - {0}
        #create decision variable vtype = Binary
        self.model.x = pyo.Var(self.model.N,self.model.N,within=pyo.Binary,initialize=0)
        #flow Constraint any customer is visited only once
        def inflow(model,i):
            return pyo.quicksum(self.model.x[i,j] for j in self.model.N if j != i) == 1
        #flow balance Constraint
        def flowBlance(model,i):
            return pyo.quicksum(self.model.x[i,j] for j in self.model.N if j !=i) ==\
                        pyo.quicksum(self.model.x[j,i] for j in self.model.N if j !=i)
        self.model.c1 = pyo.Constraint(self.model.C,rule=inflow)
        self.model.c2 = pyo.Constraint(self.model.C,rule=flowBlance)
        #Constraints of Depot
        self.model.c3 = pyo.Constraint(expr=pyo.quicksum(self.model.x[0,j] for j in self.model.C)>=1)
        self.model.c4 = pyo.Constraint(expr=pyo.quicksum(self.model.x[0,j] for j in self.model.C)<=\
                                            self.instance.vehicles)
        self.model.c5 = pyo.Constraint(expr=pyo.quicksum(self.model.x[0,j] for j in self.model.C)==\
                                            pyo.quicksum(self.model.x[j,0] for j in self.model.C))

        #Constraints Capacity
        self.model.q = pyo.Var(self.model.N,bounds=[0,self.instance.capacity],\
                               within=pyo.NonNegativeReals)
        M = self.instance.capacity
        self.model.c6 = pyo.ConstraintList()
        for i in self.model.N:
            for j in self.model.C:
                if i != j:
                    expr = self.model.q[i] + self.instance.demand[j] - (1-self.model.x[i,j])*M
                    self.model.c6.add(expr=expr<=self.model.q[j])

        #self.model.c6.deactivate()
        #self.model.c6.activate()
        #Constraints Time Windows
        self.model.t = pyo.Var(self.model.N,bounds=[0,self.instance.timeWindows[0][1]],\
                               within=pyo.NonNegativeReals)
        T = self.instance.timeWindows[0][1]
        self.model.c7 = pyo.ConstraintList()
        for i in self.model.C:
            for j in self.model.N:
                if i != j:
                    expr = self.model.t[i] + self.instance.dismatrix[i,j] + \
                           self.instance.serviceTime[i]- (1-self.model.x[i,j])*T
                    self.model.c7.add(expr=expr<=self.model.t[j])
        def twclose(model,i):
            return self.model.t[i] <= self.instance.timeWindows[i][1]
        def twopen(model,i):
            return self.model.t[i] >= self.instance.timeWindows[i][0]
        self.model.c8 = pyo.Constraint(self.model.N,rule = twopen)
        self.model.c9 = pyo.Constraint(self.model.N,rule = twclose)
        #Objective
        def objRule(model):
            return pyo.quicksum(self.instance.dismatrix[i,j] * self.model.x[i,j] \
                                for i in self.model.N for j in self.model.N)
        self.model.obj = pyo.Objective(rule=objRule,sense=pyo.minimize)

        #self.model.pprint()

    def solve(self,solver):
        opt = pyo.SolverFactory(solver)
        opt.solve(self.model)
        arcs = []
        try:
            print("最优解是：",end="")
            print("{:.2f}".format(pyo.value(self.model.obj)))
            for idx in self.model.x:
                if pyo.value(self.model.x[idx]) > 0.5:
                    arcs.append(idx)

            print("已成功获取路径信息")
            return arcs
        except:
            print("模型不可行")



    def draw(self,arcs):
        #print(arcs)
        print("开始绘制路径")
        for arc in arcs:
            x = [self.instance.coords[arc[0]][0],self.instance.coords[arc[1]][0]]
            y = [self.instance.coords[arc[0]][1],self.instance.coords[arc[1]][1]]
            plt.plot(x,y,color = "black")
        plt.show()
instance = Instance.create_instance("c101.txt")
vrptw = VRPTW(instance)
arcs = vrptw.solve(solver="cplex") #cbc gurobi
vrptw.draw(arcs=arcs)