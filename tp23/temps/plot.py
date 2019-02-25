import matplotlib.pyplot as plt
import csv

def load():
    with open('data.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        donnees = []
        for row in reader:
            #print(donnees)
            donnees.append(row)
    return donnees

def getlines(data,commande):
    res = []
    for d in data:
        if d['commande'] == str(commande):
            res.append(d)
    return res

def plot():
    data = load()
    #print(len(data))
    for commande in range(0,7):
        local = getlines(data,commande)
        plt.plot(range(0,4),[float(run["temps"]) for run in local],label="Commande #"+str(commande))
    plt.legend(range(0,7))
    plt.show()

def plot_courbes():
    plt.plot([1,2,3,4,5,6,8,10,12],[1,0.94,0.97,0.98,0.98,0.95,0.93,0.91,0.83])
    #plt.legend(range(1,13))
    plt.show()

plot_courbes()