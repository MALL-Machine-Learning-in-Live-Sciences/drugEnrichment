import decoupler as dc

net = dc.get_collectri(organism='Human')
net.to_csv('data/collectri.csv')


