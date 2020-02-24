#ifndef CRTest_PmtHit_h
#define CRTest_PmtHit_h

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

class PmtHit : public G4VHit
{
public:
	PmtHit(G4int id);
	PmtHit(const PmtHit &right);
	virtual ~PmtHit();

	const PmtHit &operator=(const PmtHit &right);
	G4bool operator==(const PmtHit &right) const;

	inline void *operator new(size_t);
	inline void operator delete(void *aHit);

	virtual void Draw();
	virtual void Print();

	G4int GetPmtID() const { return fId; }
	//void SetPmtID(G4int id){fPmtID = id;};

private:
	G4int fId;

public:
};
// HitCollection
//typedef G4THitsCollection<PmtHit> PmtHC;
using PmtHC = G4THitsCollection<PmtHit>;

extern G4ThreadLocal G4Allocator<PmtHit> *PmtHitAllocator;

inline void *PmtHit::operator new(size_t)
{
	if (!PmtHitAllocator)
	{
		PmtHitAllocator = new G4Allocator<PmtHit>;
	}
	return (void *)PmtHitAllocator->MallocSingle();
}

inline void PmtHit::operator delete(void *aHit)
{
	PmtHitAllocator->FreeSingle((PmtHit *)aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*CRTest_PmtHit_h*/