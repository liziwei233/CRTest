#include "PmtHit.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<PmtHit>* PmtHitAllocator;

PmtHit::PmtHit(G4int id)
	: fId(id)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PmtHit::~PmtHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PmtHit::PmtHit(const PmtHit &right)
: G4VHit(),
  fId(right.fId)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const PmtHit& PmtHit::operator=(const PmtHit &right)
{
  fId = right.fId;
  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool PmtHit::operator==(const PmtHit &/*right*/) const
{
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PmtHit::Draw()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PmtHit::Print()
{
  G4cout << "  PmtHit[" << fId << "] " << G4endl;	
}