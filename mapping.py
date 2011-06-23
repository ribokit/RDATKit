import random
import settings

class MappingData:
    def __init__(self, data=[], seqpos=[], type=''):
	if seqpos:
	    self._data = [None]*(max(seqpos) + 1)
	    self.seqpos = seqpos
	    for i, pos in enumerate(seqpos):
		self._data[pos] = data[i]
	else:
	    self._data = data
	    self.seqpos = range(len(data))
	self.type = type
    
    def load(self, shapefile):
        self.seqpos = []
	mdata = []
        for line in shapefile.readlines():
	    fields = line.strip().split(' ')
	    self.seqpos.append(int(fields[0])-1)
	    mdata.append(float(fields[1]))
	self._data = [None]*(max(self.seqpos) + 1)
	for i, dat in enumerate(mdata):
	    self._data[self.seqpos[i]] = dat
	    

    def __getitem__(self, k):
        return self._data[k]

    def __str__(self):
        s = ''
        for pos in self.seqpos:
	    if self._data[pos] is not None:
		s += '%d %d\n' % (pos + 1, float(self._data[pos]))
        return s
    
    def sample(self, numsamples, replacement=False):
        if replacement:
	    nseqpos = [0]*numsamples
	    for i in range(numsamples):
		idx = random.choice(self.seqpos)
		nseqpos[i] = idx
	else:
	    nseqpos = random.sample(self.seqpos, numsamples)
	ndata = [None]*len(nseqpos)
	for i, pos in enumerate(nseqpos):
	    ndata[i] = self._data[pos]
	return MappingData(data=ndata, seqpos=nseqpos, type=self.type)
	     
