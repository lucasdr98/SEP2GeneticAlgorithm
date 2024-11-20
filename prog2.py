import numpy as np
from geneticalgorithm import geneticalgorithm as ga

# Classe Linha
class Linha():
    def __init__(self, de, para, res, reat, SUCsh, TAP, DEF, status, capacidade):
        self.de = de
        self.para = para
        self.res = res
        self.reat = reat
        self.SUCsh = SUCsh
        self.TAP = TAP
        self.DEF = DEF
        self.status = status
        self.capacidade = capacidade
        self.fluxoP = 0
        self.fluxoPmk = 0
        self.fluxoQ = 0
        self.fluxoQmk = 0
        self.perdarP = 0
        self.perdarQ = 0
    def getY(self):
        Y = 1/(self.res + self.reat)
        return Y

# Classe Barra
class Barra():
    def __init__(self, numero, tipo, tensao, angulo, linhas, PD, QD, Bsh, PG, custo,CGMin, CGMax):
        self.numero = numero
        self.tipo = tipo
        self.tensao = tensao
        self.angulo = angulo
        self.linhas = linhas
        self.PD = PD
        self.QD = QD
        self.Bsh = Bsh
        self.PG = PG
        self.custo = custo
        self.CGMax = CGMax
        self.CGMin = CGMin
        self.Pcalc = 0
        self.Qcalc = 0
    
    def pEsp(self):
        return self.PG - self.PD
     

# Função para ler o arquivo e criar objetos das classes
def lerArquivo(nomeArquivo):
    barras = []
    lts = []
    EmDBAR = False
    EmDlinha = False
    stop = 0
    
    arquivo = open(nomeArquivo, 'r')
    linhas = arquivo.readlines()
    arquivo.close()
    linhas = [linha.replace('\n', '') for linha in linhas]
    linhas = [linha.split(';') for linha in linhas]

    for linha in linhas:
        #Filtragem das linhas

        #Verifica se esta na linha de dados das barras
        if 'DBAR' in linha:
            EmDlinha = False
            EmDBAR = True
            continue
        #Verifica se esta na linha de dados das linhas
        if 'DCIR' in linha:
            EmDlinha = True
            EmDBAR = False
            continue
        #ignora a linha de cabeçalho
        if linha[0].startswith('BARRA') or linha[0].startswith(' BDE'):
            continue
        #ignora as linhas com x
        if linha[0].startswith('x'):
            continue
        #ignora as linhas com #
        if linha[0].startswith('#'):
            stop += 1
            if stop == 2:
                return barras, lts
            else:
                continue
        #ignora as linhas vazias
        if linha == ['']:
            continue

        linha = linha[0].split()
        #Cria os objetos
        if EmDBAR:
            barra_numero = int(linha[0])
            barra_PD = float(linha[1])
            barra_QD = float(linha[2])
            barra_Bsh = complex(0,float(linha[3]))
            barra_Tipo = linha[4]
            barra_Tensao = float(linha[5])
            barra_Angulo = np.deg2rad(float(linha[6]))
            barra_PG = float(linha[7])
            barra_Custo = float(linha[8])
            barra_CGMin = float(linha[9])
            barra_CGMax = float(linha[10])
            barra = Barra(barra_numero, barra_Tipo, barra_Tensao, barra_Angulo, [], barra_PD, barra_QD, barra_Bsh, barra_PG, barra_Custo,barra_CGMin, barra_CGMax)
            barras.append(barra)

        if EmDlinha:
            lt_de = int(linha[0])
            lt_para = int(linha[1])
            lt_res = float(linha[3])
            lt_reat = complex(0,float(linha[4]))
            lt_SUCsh = complex(0,float(linha[5]))
            lt_TAP = float(linha[6])
            lt_DEF = float(linha[7])
            lt_status = linha[8]
            lt_capacidade = float(linha[9])
            lt = Linha(lt_de, lt_para, lt_res, lt_reat, lt_SUCsh, lt_TAP, lt_DEF, lt_status, lt_capacidade)
            lts.append(lt)
            #adiciona os objetos das linhas ao objeto das barras
            for barra in barras:
                if (barra.numero == lt_de) or (barra.numero == lt_para):
                    barra.linhas.append(lt)
# Função para gerar o Ybus
def geraYBUS(linhas, barras):
    Ybus = np.zeros((len(barras), len(barras)), dtype=complex)
    
    for linha in linhas:
        if linha.status == 'L':
            k = linha.de - 1
            m = linha.para - 1

            Ybus[k][k] += linha.getY() * (linha.TAP**2) + (linha.SUCsh/2)
            Ybus[m][m] += linha.getY() + (linha.SUCsh/2)
            Ybus[k][m] += -linha.getY()*linha.TAP*(np.cos(-linha.DEF) + 1j*np.sin(-linha.DEF))
            Ybus[m][k] += -linha.getY()*linha.TAP*(np.cos(linha.DEF) + 1j*np.sin(linha.DEF))

    for barra in barras:
        Ybus[barra.numero - 1][barra.numero - 1] += barra.Bsh

    return Ybus
# Função para calcular os deltas
def calculaDeltas(barras, Ybus):
    G = np.real(Ybus)
    B = np.imag(Ybus)
    DeltaP = []
    DeltaQ = []
    DeltaPQ = []
    for barraK in barras:
        P = 0
        Q = 0
        for barraM in barras:
            k = barraK.numero - 1
            m = barraM.numero - 1
            Theta = barraK.angulo - barraM.angulo
            
            Q += barraM.tensao*(G[k][m]*np.sin(Theta) - B[k][m]*np.cos(Theta))
            P += barraM.tensao*(G[k][m]*np.cos(Theta) + B[k][m]*np.sin(Theta))
            
            
        P = P*barraK.tensao
        Q = Q*barraK.tensao

        barraK.Qcalc = Q
        barraK.Pcalc = P
        
        if (barraK.tipo != 'PQ'):
            DeltaQ.append(0)
        else:
            DeltaQ.append(-barraK.QD - Q)

        if (barraK.tipo == 'SW'):
            DeltaP.append(0)
        else:
            DeltaP.append(barraK.pEsp() - P)
        
        
        
    DeltaPQ = np.concatenate((DeltaP, DeltaQ))
    return DeltaPQ
# Função para calcular o Jacobiano
def Jacobiana(barras, Ybus):
    G = np.real(Ybus)
    B = np.imag(Ybus)

    H = np.zeros((len(barras), len(barras)))
    N = np.zeros((len(barras), len(barras)))
    M = np.zeros((len(barras), len(barras)))
    L = np.zeros((len(barras), len(barras)))
    

    
    for barraK in barras:
        somaH = 0
        SomaN = 0
        SomaM = 0
        SomaL = 0

        for barraM in barras:
            Theta = barraK.angulo - barraM.angulo
            k = barraK.numero - 1
            m = barraM.numero - 1

            #Calculando H
            H[k][m] = barraK.tensao*barraM.tensao*(G[k][m]*np.sin(Theta) - B[k][m]*np.cos(Theta))
            somaH += barraM.tensao*(G[k][m]*np.sin(Theta)-B[k][m]*np.cos(Theta))

            #Calculando N
            N[k][m] = barraK.tensao*(G[k][m]*np.cos(Theta)+B[k][m]*np.sin(Theta))
            SomaN += barraM.tensao*(G[k][m]*np.cos(Theta)+B[k][m]*np.sin(Theta))

            #Calculando M
            M[k][m] = -barraK.tensao*barraM.tensao*(G[k][m]*np.cos(Theta)+B[k][m]*np.sin(Theta))
            SomaM += barraM.tensao*(G[k][m]*np.cos(Theta)+B[k][m]*np.sin(Theta))

            #Calculando L
            L[k][m] = barraK.tensao*(G[k][m]*np.sin(Theta) - B[k][m]*np.cos(Theta))
            SomaL += barraM.tensao*(G[k][m]*np.sin(Theta) - B[k][m]*np.cos(Theta))

        H[k][k] = -np.power(barraK.tensao, 2) * B[k][k] - barraK.tensao * somaH
        N[k][k] = barraK.tensao * G[k][k] + SomaN
        M[k][k] = -np.power(barraK.tensao, 2) * G[k][k] + barraK.tensao * SomaM
        L[k][k] = -barraK.tensao * B[k][k] + SomaL
    
    #Tratando as matrizes
    for barraK in barras:
        k = barraK.numero - 1
        for barraM in barras:
            m = barraM.numero - 1

            if (barraK.tipo != 'PQ'):
                L[k][m] = 0
                L[m][k] = 0
                N[m][k] = 0
                M[k][m] = 0
            if (barraK.tipo == 'SW'):
                H[k][m] = 0
                H[m][k] = 0
                M[m][k] = 0
                N[k][m] = 0
        
        if (barraK.tipo != 'PQ'):
            L[k][k] = float("inf")

        if (barraK.tipo == 'SW'):
            H[k][k] = float("inf")
    J = np.block([[H, N], [M, L]])
    return J
# Função para aplicar o método de Newton-Raphson             
def newtonRaphson(barras, Ybus):
    tolerancia = 0.0001
    erro = 1000
    DeltaPQ = calculaDeltas(barras, Ybus) # Calculando os deltas iniciais
    ite = 0 #contador de iterações
    while erro > tolerancia:
        J = Jacobiana(barras, Ybus) # Calculando o Jacobiano
        
        Thetas = []
        Vs = []
        #Montando x
        for barra in barras:
            Thetas.append(barra.angulo)
            Vs.append(barra.tensao)
        x = np.concatenate((Thetas, Vs))

        #Calculando Deltax
        Deltax = np.matmul(np.linalg.inv(J), DeltaPQ)
        x = x + Deltax

        #Substituindo valores nas barras
        newThetaV = np.array_split(x, 2)
        for barra in barras:
            barra.angulo = newThetaV[0][barra.numero - 1]
            barra.tensao = newThetaV[1][barra.numero - 1]

        DeltaPQ = calculaDeltas(barras, Ybus) # Calculando os deltas
        erro = np.max(np.abs(DeltaPQ)) # Calculando o erro
        ite += 1 # Incrementando o contador

        if ite > 20:
            return - 1
    
    return barras
# Função para calcular os fluxos
def calculaFluxos(barras, linhas):
    Thetas = []
    Vs = []
    #Pegando valores de tensão e angulo
    for barra in barras:
            Thetas.append(barra.angulo)
            Vs.append(barra.tensao)

    for linha in linhas:
        k = linha.de - 1
        m = linha.para - 1
        y = linha.getY()
        g = np.real(y)
        b = np.imag(y)
        bsh = (np.imag(linha.SUCsh)/2)
        Thetakm = Thetas[k] - Thetas[m]
        Thetamk = Thetas[m] - Thetas[k]

        Pkm = ((linha.TAP*Vs[k])**2)*g - (linha.TAP*Vs[k])*Vs[m]*g*np.cos(Thetakm+linha.DEF) - (linha.TAP*Vs[k])*Vs[m]*b*np.sin(Thetakm+linha.DEF)
        Pmk = (Vs[m]**2)*g - (linha.TAP*Vs[k]*Vs[m]*g*np.cos(Thetamk - linha.DEF)) - (linha.TAP*Vs[k]*Vs[m]*b*np.sin(Thetamk - linha.DEF))
        
        Qkm = -((linha.TAP*Vs[k])**2)*(b + bsh) + (linha.TAP*Vs[k])*Vs[m]*b*np.cos(Thetakm + linha.DEF) - (linha.TAP*Vs[k])*Vs[m]*g*np.sin(Thetakm + linha.DEF)
        Qmk = -(Vs[m]**2)*(b + bsh) + (linha.TAP*Vs[k]*Vs[m]*b*np.cos(Thetamk - linha.DEF)) - (linha.TAP*Vs[k]*Vs[m]*g*np.sin(Thetamk - linha.DEF))
        
        PerdaP = Pkm + Pmk
        PerdaQ = Qkm + Qmk

        linha.perdaP = PerdaP
        linha.perdaQ = PerdaQ

        linha.fluxoP = Pkm
        linha.fluxoPmk = Pmk
        linha.fluxoQ = Qkm
        linha.fluxoQmk = Qmk

    return linhas
# Função para salvar os resultados
def salvaResultados(barras, linhas):
    with open('resultadosEC2.txt', 'w', encoding="utf-8") as arquivo:
        arquivo.write('Resultados barras\n')
        arquivo.write('Barra   Tensão   Angulo      PI        PG        QI        QG' + '\n')
        for barra in barras:
            
            arquivo.write(f'{barra.numero}      {barra.tensao: .4f}  {np.rad2deg(barra.angulo): .4f}   {barra.Pcalc: .4f}   {barra.PD+barra.Pcalc: .4f}   {barra.Qcalc+(np.imag(barra.Bsh)*barra.tensao**2): .4f}   {barra.QD+barra.Qcalc: .4f}' + '\n')

        arquivo.write('\n')
        arquivo.write('\n')

        arquivo.write('Resultados linhas\n')
        arquivo.write('de   para  fluxoPkm fluxoPmk   fluxoQkm  fluxoQmk   fluxoS   perdaP    perdaQ' + '\n')
        for linha in linhas:
            arquivo.write(f'{linha.de}     {linha.para}   {linha.fluxoP: .4f}   {linha.fluxoPmk: .4f}   {linha.fluxoQ: .4f}    {linha.fluxoQmk: .4f}   {np.sqrt(linha.fluxoP**2 + linha.fluxoQ**2): .4f}  {linha.perdaP: .4f}   {linha.perdaQ: .4f}' + '\n')

def checaLimites(VAR):
    barras, linhas = lerArquivo('dados_sistema12B_EC2_CasoBase.txt') # lendo o arquivo e criando objetos
    # inserindo valores nas barras e linhas
    barras[0].PG = VAR[0]
    barras[3].PG = VAR[1]
    barras[8].PG = VAR[2]
    linhas[7].TAP = VAR[3]
    linhas[8].TAP = VAR[4]
    linhas[7].DEF = np.deg2rad(VAR[5])
    linhas[8].DEF = np.deg2rad(VAR[6])
    barras[0].tensao = VAR[7]
    barras[2].tensao = VAR[8]
    barras[3].tensao = VAR[9]
    barras[8].tensao = VAR[10]
    
    Ybus = geraYBUS(linhas, barras) # Gerando a Ybus
    barras = newtonRaphson(barras, Ybus)
    if barras == -1:  # Checa se o NR convergiu
        return 100000 

    linhas = calculaFluxos(barras, linhas)

    SumDelta = 0

    # Checando limites das barras
    for barra in barras:

        # Checando limite de tensão
        if barra.tipo == 'PQ':  
            if barra.numero <= 6: #Filtrando as barras de 500kV
                if barra.tensao > 1.05 or barra.tensao < 1:
                    SumDelta += 10*np.abs(barra.tensao - 1)
                else:
                    SumDelta += np.abs(barra.tensao - 1)  # Somatorio das diferenças de tensão para 1 p.u.
                
            else:
                if barra.tensao > 1.05 or barra.tensao < 0.96:
                    SumDelta += 10*np.abs(barra.tensao - 1)
                else:
                    SumDelta += np.abs(barra.tensao - 1)  # Somatorio das diferenças de tensão para 1 p.u.

        # Checando limite de potencia
        if (barra.tipo == 'SW'):
            if (barra.PD+barra.Pcalc) > 0.9*barra.CGMax or (barra.PD+barra.Pcalc) < barra.CGMin:
                SumDelta += 10 * np.abs((barra.PD+barra.Pcalc) - barra.CGMax)
            else:
                SumDelta += (barra.PD+barra.Pcalc - barra.CGMax)
    

    # Checando limites das linhas
    for linha in linhas:
        Skm = np.sqrt(linha.fluxoP**2 + linha.fluxoQ**2)
        Smk = np.sqrt(linha.fluxoPmk**2 + linha.fluxoQmk**2)

        if Skm > 0.9*linha.capacidade:
            SumDelta += 10*(Skm - linha.capacidade)
        else:
            SumDelta += (Skm - linha.capacidade)  # Somatorio das diferenças de capacidade

        if Smk > 0.9*linha.capacidade:
            SumDelta += 10*(Smk - linha.capacidade)
        else:
            SumDelta += (Smk - linha.capacidade)  # Somatorio das diferenças de capacidade
        
        if Skm < 0.01 or Smk < 0.01: # Verifica se alguma linha esta deixando de ser usada
            SumDelta += 10
        
        

    return SumDelta


# Limites para G1, G4, G9, TAP5, TAP6, DEF5, DEF6, V1, V3, V4, V9
bounds = np.array([(0.2, 1.5), (0.2, 1.4), (0.2, 1.1), (0.2, 1.0), (0.2, 1.0), (-30, 30), (-30, 30), (1, 1.05), (1, 1.05), (1, 1.05), (0.95, 1.05)])

ga_instance = ga(function=checaLimites, 
                 dimension=11, 
                variable_type='real', 
                variable_boundaries=bounds, 
                algorithm_parameters={'max_num_iteration': None,
                                        'population_size': 300,
                                        'mutation_probability': 0.05,
                                        'elit_ratio': 0.2,
                                        'crossover_probability': 0.8,
                                        'parents_portion': 0.3,
                                        'crossover_type': 'uniform',
                                        'max_iteration_without_improv': 30
                                        },
                function_timeout = 60,
                convergence_curve=True)

ga_instance.run()
solution = ga_instance.best_variable

print(solution)

barras, linhas = lerArquivo('dados_sistema12B_EC3_CasoBase.txt') # lendo o arquivo e criando objetos

     
# # inserindo valores nas barras e linhas
barras[0].PG = solution[0]
barras[3].PG = solution[1]
barras[8].PG = solution[2]
linhas[7].TAP = solution[3]
linhas[8].TAP = solution[4]
linhas[7].DEF = np.deg2rad(solution[5])
linhas[8].DEF = np.deg2rad(solution[6])
barras[0].tensao = solution[7]
barras[2].tensao = solution[8]
barras[3].tensao = solution[9]
barras[8].tensao = solution[10]
    
Ybus = geraYBUS(linhas, barras) # Gerando a Ybus
barras = newtonRaphson(barras, Ybus)
print(barras)
linhas = calculaFluxos(barras, linhas)

salvaResultados(barras, linhas)
