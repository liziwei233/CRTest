
#include "Inc/CRInfo.h"
#include "Inc/hit_data_dec_total.hh"
using namespace std;

//class CRTracker;
//class CRT0;

//gSystem->Load("hit_data_dec_total.cc");
TString path = "/mnt/d/Experiment/DATA/cosmicray/20210207";

/* 3.1 data
int Eabled[FTOFChNum] = {1, 1, 1, 1,
                         1, 0, 0, 1,
                         0, 0, 0, 1,
                         1, 1, 1, 1};
int FTOFELEchid[] = {0, 1, 2, 3,
                     4, 5, 6, 7,
                     8, 9, 10, 11,
                     12, 13, 14, 15};
int FTOFDetchid[] = {0, 1, 2, 3,
                     8, -1, -1, 7,
                     -1, -1, -1, 11,
                     12, 13, 14, 15};
int T0ELEchid[] = {8, 9, 10};
int T0Detchid[] = {0, 1, 2};
*/
/* 3.3-3.5 data
int Eabled[FTOFChNum] = {1, 1, 1, 1,
                         1, 1, 1, 1,
                         0, 0, 0, 1,
                         1, 1, 1, 1};
int FTOFELEchid[] = {0, 1, 2, 3,
                     4, 5, 6, 7,
                     8, 9, 10, 11,
                     12, 13, 14, 15};
int FTOFDetchid[] = {0, 1, 2, 3,
                     8, 9, 10, 7,
                     -1, -1, -1, 11,
                     12, 13, 14, 15};
int T0ELEchid[] = {8, 9, 10};
int T0Detchid[] = {2, 1, 0};
*/

const int NumMM = 4;
void SetPath(TString YourPath)
{
    path = YourPath;
};
// **
// * TODO *

void ReadTrackAGTData2Root(vector<TString> datList, TString fRawName, int force)
{
    if (datList.size() == 0)
    {
        cout << "Error!Dat size=0!" << endl;
        return;
    }
    ofstream op = LogFile(fRawName, "ReadTrackAGTData2Root");
    int length = sizeof(unsigned short);
    unsigned short memblock;
    int boardstart = 0;
    int chipstart = 0;
    int type;
    std::string filename;
    bool headerFlag;
    bool trailerFlag;
    double counting;
    unsigned short triggerType;
    int badNumber = 0;
    int insidebadNum = 0;
    int preevent;
    int event;
    int carry;
    int breaknum;
    int carrynum;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> wave;
    vector<int> pos;
    double wave_max;
    double charge;
    double ped;
    double ped_rms;

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fRawName << " is exist!" << endl;
            return;
        }
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
    {
        cout << "ERROR! " << fRawName << " cant open!" << endl;
        return;
    }
    TTree *fTree = new TTree("tree", "tree");
    fTree->Branch("event", &event);
    fTree->Branch("board", &board);
    fTree->Branch("chip", &chip);
    fTree->Branch("channel", &channel);
    fTree->Branch("wave", &wave);
    fTree->Branch("pos", &pos);

    TTree *fTree2 = new TTree("tree2", "tree2");
    fTree2->Branch("type", &type);
    type = AGETType;
    fTree2->Fill();
    fTree2->Write();

    TTree *fTree3 = new TTree("tree3", "tree3");
    fTree3->Branch("filename", &filename);
    vector<int> bdList;
    TH1I *hid = new TH1I("hid", "hid", 100e4, 1, 100e4 + 1);
    preevent = 0;
    carry = 0;
    breaknum = 0;
    carrynum = 0;
    for (int i = 0; i < (int)datList.size(); i++)
    {
        insidebadNum = 0;
        fstream InDat;
        InDat.open(datList[i], ios::in | ios::binary);
        cout << "--> Reading a new dat file-" << i << " : " << datList[i] << endl;
        filename = GetFilename(datList[i]);
        fTree3->Fill();

        while (!InDat.eof())
        {
            InDat.read((char *)(&memblock), length);
            memblock = AGETExchange(memblock);
            if (memblock == 0xeeee)
            {
                headerFlag = 1;
                counting = 0;
                event = 0; // reset
                           //cout << "header" << endl;
            }              // data packet header
            //else if (memblock&0xe == 0xe && headerFlag == 1 && counting == 0)	// board number, chip number, trigger type
            else if (counting == 1) // board number, chip number, trigger type
            {
                //board = (((memblock&0x0f00)>>8 - boardstart)); //cout << "board number = " << board << endl;
                board = (((memblock & 0x0f00) >> 8) - boardstart); //cout << "board number = " << board << endl;
                //if (board == 7) board = 5;
                //if (board == 9) board = 1;
                if (bdList.size() < 4)
                {
                    int flag = -1;
                    for (size_t ii = 0; ii < bdList.size(); ii++)
                        if (board == bdList[ii])
                            flag = 1;
                    if (flag == -1)
                        bdList.push_back(board);
                }

                chip = (((memblock & 0x00f0) >> 4) - chipstart); //cout << "chip number = " << chip << endl;
                triggerType = (memblock & 0x000f);               //cout << "trigger type = " << triggerType << endl;
                                                                 //cout << "fixed number" << endl;// counting = 1;
            }
            else if (counting > 1 && counting < 6)
            {
            }
            else if (counting == 6)
            {
                channel = memblock; /*cout << "channel number = " << channel << endl;*/
            }
            else if (counting > 6 && counting < 10)
            {
                if (counting == 7)
                    event += memblock * (2 ^ 32);
                else if (counting == 8)
                    event += memblock * (2 ^ 16);
                else if (counting == 9)
                {
                    event += memblock; /*cout << "trigger Number = " << event << endl;*/
                }
            }
            //else if (counting > 9 && counting < 521)
            else if (memblock == 0xffff)
            {
                trailerFlag = 1;
                if (wave.size() == 512)
                {

                    wave_max = TMath::MaxElement(300, &wave.at(100));
                    carry = IsCarry(event, preevent);
                    if (carry)
                        continue;
                    /*
                    if (carry == 1e4)
                        breaknum++;
                    if (carry == 65535)
                        carrynum++;
                    if (0)
                    {
                        op << "--> Reading a new dat file-" << i << " : " << datList[i] << endl;
                        op << "event, preevent, carry: " << event << "," << preevent << "," << carry << endl;
                        op << "breaknum, carrynum: " << breaknum << ", " << carrynum << endl;
                        op << "calculated event:" << breaknum * 1e4 + carrynum * 65535 + event << endl;
                        op << endl;
                    }
                    */
                    preevent = event;
                    event = breaknum * 1e4 + carrynum * 65535 + event;

                    fTree->Fill();
                    hid->Fill(event);
                }
                else if (wave.size() != 0)
                { /*cout << "wave.size() = " << wave.size() << endl;*/
                    insidebadNum++;
                    badNumber++;
                    event = 0;
                    //preevent = event;
                }
                wave.clear();
                pos.clear();
                //cout << "clear ..." << endl;
            }
            else
            {
                double buffer = memblock & 0x0fff;
                wave.push_back(buffer);
                pos.push_back(wave.size());
                //cout << counting << "\t" << buffer << "\t" << wave.size() << endl;
            }

            counting++;
        }
        if (insidebadNum > 0)
        {
            cout << " [ Warning! ] Dat file-" << i << " : " << datList[i] << "have BAD EVENTS=" << insidebadNum << "!" << endl;
        }
    }
    fTree3->Write();

    cout << "--> Total entries = " << fTree->GetEntries() << endl;
    cout << "--> Bad Events = " << badNumber << endl;
    cout << "--> Board list: ";
    for (size_t i = 0; i < bdList.size(); i++)
        cout << bdList[i] << " ";
    cout << endl;
    //hid->Draw();
    int Record = 0;
    int tristart = 0;
    int triend = 0;
    bool flag = 0;

    for (int i = hid->FindFirstBinAbove(63); i < hid->FindLastBinAbove(63) + 1; i++)
    {

        if (hid->GetBinContent(i) >= 63 && hid->GetBinContent(i - 1) >= 63 && hid->GetBinContent(i - 2) >= 63)
        {
            if (!flag)
            {
                tristart = i - 2;
                flag = 1;
            }
            else
                triend = i;
        }
    }
    for (int i = tristart; i < triend + 1; i++)
    {
        if (hid->GetBinContent(i) >= 63)
            Record++;
    }

    cout << "--> The start id: " << tristart << ", end id: " << triend << ", Total Triggers = " << triend - tristart + 1 << endl;
    cout << "--> Recorded Events = " << Record << endl;
    fTree->Write();
    //fTree->Print();
    fFile->Flush();
    fFile->Close();
    //return;

    cout << "--> Track-AGET Raw data file has been convert to root-file: " << fRawName << ".\n"
         << endl;
};

void AnalysisAGET(TString fRawName, TString fAnaName, int force)
{

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fAnaName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fAnaName << " is exist!" << endl;
            return;
        }
    }

    TFile *anaFile = new TFile(fAnaName, "recreate");
    if (!anaFile->IsOpen())
    {
        cout << "ERROR! " << fAnaName << " cant open!" << endl;
        return;
    }

    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;

    TTree *htree[4];
    AGETHit dethit;
    for (int i = 0; i < 4; i++)
    {

        htree[i] = new TTree(Form("tree%d", i), Form("tree%d", i));
        htree[i]->Branch("dethit", &dethit);
    }

    double factor = 2e-7 / 50 / 1e3 / 1e-9; // factor that converts to nC.
    double deltat = 200;                    // Unit:ns, 5M sample rate
    TFile *rawFile = new TFile(fRawName, "read");
    cout << "--> Analyzing waveform files from " << fRawName << endl;
    TH1D *hNoise = new TH1D("hNoise", "hNoise", 1000, 0, 1000);
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);
    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;
        if (ii % 1000 == 0)
            cout << " --> Entry$= " << ii << endl;
        tree->GetEntry(ii);
        hNoise->Reset();
        dethit.event = 0;
        dethit.q = 0;
        dethit.t = 0;
        dethit.wave_max = 0;
        dethit.risetime = 0;
        dethit.ped = 0;
        dethit.ped_rms = 0;
        dethit.chn[0] = -999;
        dethit.chn[1] = -999;
        for (int j = 0; j < (int)wave->size(); j++)
        {
            if ((j > 5 && j < 150) || (j > 320 && j < 500))
            {

                hNoise->Fill(wave->at(j));
            }
            /*
            if (j >= 150 && j < 320)
            {
                dethit.q += wave->at(j) * factor;
            }
            */
        }
        dethit.ped = hNoise->GetMean();
        dethit.ped_rms = hNoise->GetRMS();
        if (hNoise->GetEntries() > 0)
        {
            int status = 1;
            if (dethit.ped - 3 * dethit.ped_rms < 50)
                status = hNoise->Fit("gaus", "RQ", "", 50, dethit.ped + 3 * dethit.ped_rms);
            else
                status = hNoise->Fit("gaus", "RQ", "", dethit.ped - 3 * dethit.ped_rms, dethit.ped + 3 * dethit.ped_rms);
            if (status == 0)
            {

                dethit.ped = hNoise->GetFunction("gaus")->GetParameter(1);
                dethit.ped_rms = hNoise->GetFunction("gaus")->GetParameter(2);
            }
            else
                cout << ">> ERROR! Fit false! <<" << endl;
        }
        //dethit.q = dethit.q - dethit.ped * (320 - 150) * factor;
        dethit.wave_max = TMath::MaxElement(320 - 150, &(wave->at(150))) - dethit.ped;
        dethit.t = TMath::LocMax(320 - 150, &(wave->at(150))) + 150;
        dethit.event = event;
        vector<double> waveform;
        waveform.clear();
        for (int i = 0; i < wave->size(); i++)
        {
            waveform.push_back(wave->at(i) - dethit.ped);
        }
        double cf = 0.1 * dethit.wave_max;
        double risetime1 = 0;
        double risetime2 = 0;
        for (int i = dethit.t; i >= 1 && i < 512; --i)
        {
            if (waveform[i] - cf > 0 && waveform[i - 1] - cf < 0)
            {
                double x1 = i - 1;
                double x2 = i;
                double y1 = waveform[i - 1];
                double y2 = waveform[i];
                risetime1 = linear_interpolation(x1, x2, y1, y2, cf);
                break;
            }
        }

        cf = 0.9 * dethit.wave_max;
        for (int i = dethit.t; i >= 1 && i < 512; --i)
        {
            if (waveform[i] - cf > 0 && waveform[i - 1] - cf < 0)
            {
                double x1 = i - 1;
                double x2 = i;
                double y1 = waveform[i - 1];
                double y2 = waveform[i];
                risetime2 = linear_interpolation(x1, x2, y1, y2, cf);
                break;
            }
        }
        dethit.risetime = (risetime2 - risetime1) * deltat;

        for (int j = 0; j < 2; j++)
        {
            if (IndexXY(board, chip) == j)
            {
                // y axis
                dethit.chn[j] = channel;
            }
        }

        //cout <<"board "<<board<<endl;
        //cout <<"chip "<<chip<<endl;
        //cout <<"IndexDet(board, chip) "<<IndexDet(board, chip)<<endl;
        htree[IndexDet(board, chip)]->Fill();
    }
    for (int i = 0; i < 4; i++)
    {
        anaFile->WriteTObject(htree[i]);
    }
    anaFile->Close();
    rawFile->Close();
    cout << "--> Analyzed results has been stored to " << fAnaName << endl
         << endl;
}

void CheckTrackerID(TString fAnaName)
{
    cout << "open your file: " << fAnaName << endl;
    //TFile *hitFile = new TFile(fHitName, "recreate");
    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[4];
    AGETHit *dethit;
    dethit = new AGETHit;
    TString fpath = GetPath(fAnaName);
    ofstream op(Form("%s/Det_eventinfo.dat", fpath.Data()));

    vector<int> idvecx[4];
    vector<int> idvecxpos[4];
    vector<int> idvecx_entrynum[4];
    vector<int> idvecy[4];
    vector<int> idvecypos[4];
    vector<int> idvecy_entrynum[4];
    TGraph *gx[4];
    TGraph *gy[4];

    for (int s = 0; s < 4; s++)
    {
        cout << "Entry Detector " << s << " ..." << endl;
        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        tree[s]->SetBranchAddress("dethit", &dethit);

        idvecx[s].clear();
        idvecy[s].clear();
        idvecxpos[s].clear();
        idvecypos[s].clear();

        int nentries = tree[s]->GetEntriesFast();
        cout << "Entries is " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {
            Long64_t ii = tree[s]->LoadTree(i);
            if (ii < 0)
                break;

            tree[s]->GetEntry(ii);
            //cout<<"event id: "<<dethit->event<<endl;
            // x axis has a hit
            //if (dethit->wave_max > 7 * dethit->ped_rms && dethit->q > 5)
            if (1)
            {

                if (dethit->chn[0] != -999)
                {
                    if (idvecx[s].size() < 1)
                    {

                        idvecx[s].push_back(dethit->event);
                        idvecxpos[s].push_back(idvecx[s].size());
                        idvecx_entrynum[s].push_back(i);
                    }
                    if (idvecx[s].back() != dethit->event)
                    {

                        idvecx[s].push_back(dethit->event);
                        idvecxpos[s].push_back(idvecx[s].size());
                        idvecx_entrynum[s].push_back(i);
                    }
                    //cout<<"x (board2) record event"<<idvecx[s].back()<<endl;
                }

                // y axis has a hit
                if (dethit->chn[1] != -999)
                {
                    if (idvecy[s].size() < 1)
                    {

                        idvecy[s].push_back(dethit->event);
                        idvecypos[s].push_back(idvecy[s].size());
                        idvecy_entrynum[s].push_back(i);
                    }
                    if (idvecy[s].back() != dethit->event)
                    {

                        idvecy[s].push_back(dethit->event);
                        idvecypos[s].push_back(idvecy[s].size());
                        idvecy_entrynum[s].push_back(i);
                    }
                    // cout<<"y (board1) record event"<<idvecy.back()<<endl;
                }
            }
        }

        int xlostcounter = CalculateLostID(idvecx[s], Form("%s/Xlostevent.dat", fpath.Data()));
        int ylostcounter = CalculateLostID(idvecy[s], Form("%s/Ylostevent.dat", fpath.Data()));
        int xweiredcounter = CalculateWeiredID(idvecx[s], Form("%s/XweiredNum.dat", fpath.Data()));
        int yweiredcounter = CalculateWeiredID(idvecy[s], Form("%s/YweiredNum.dat", fpath.Data()));

        /*
        for (int i = 1; i < idvecx.size() - 1; i++)
        {
            if (idvecx.at(i) - idvecx.at(i - 1) != 1)

                cout << "x (board2): " << i << ",\t" << idvecx.at(i - 1) << ",\t" << idvecx.at(i) << ",\t" << idvecx.at(i + 1) << endl;
        }
        for (int i = 1; i < idvecy.size() - 1; i++)
        {
            if (idvecy.at(i) - idvecy.at(i - 1) != 1)
                cout << "y (board1): " << idvecy.at(i - 1) << ",\t" << idvecy.at(i) << ",\t" << idvecy.at(i + 1) << endl;
        }
*/
        gx[s] = new TGraph(idvecx[s].size(), &idvecxpos[s][0], &idvecx[s][0]);
        gy[s] = new TGraph(idvecy[s].size(), &idvecypos[s][0], &idvecy[s][0]);
        cout << "-->  The tacker " << s << " event id has been completed: " << endl
             << fAnaName << endl
             << "    |-->  The total numbers of Tracker event (x): " << idvecx[s].size() << endl
             << "    |-->  The lost  numbers of Tracker event (x): " << xlostcounter << ", eff: " << (double)xlostcounter / idvecx[s].size() << endl
             << "    |-->  The weired numbers of Tracker event (x): " << xweiredcounter << ", eff: " << (double)xweiredcounter / idvecx[s].size() << endl

             << "    |-->  The total numbers of Tracker event (y): " << idvecy[s].size() << endl
             << "    |-->  The lost  numbers of Tracker event (x): " << ylostcounter << ", eff: " << (double)ylostcounter / idvecx[s].size() << endl
             << "    |-->  The weired numbers of Tracker event (x): " << yweiredcounter << ", eff: " << (double)yweiredcounter / idvecx[s].size() << endl
             << endl;
        op << "-->  The tacker " << s << " event id has been completed: " << endl
           << fAnaName << endl
           << "    |-->  The total numbers of Tracker event (x): " << idvecx[s].size() << endl
           << "    |-->  The lost  numbers of Tracker event (x): " << xlostcounter << ", eff: " << (double)xlostcounter / idvecx[s].size() << endl
           << "    |-->  The weired numbers of Tracker event (x): " << xweiredcounter << ", eff: " << (double)xweiredcounter / idvecx[s].size() << endl

           << "    |-->  The total numbers of Tracker event (y): " << idvecy[s].size() << endl
           << "    |-->  The lost  numbers of Tracker event (x): " << ylostcounter << ", eff: " << (double)ylostcounter / idvecx[s].size() << endl
           << "    |-->  The weired numbers of Tracker event (x): " << yweiredcounter << ", eff: " << (double)yweiredcounter / idvecx[s].size() << endl
           << endl;
    }
    TLatex *ll;
    TCanvas *c = new TCanvas("c5", "c5", 1600, 1000);
    c->Divide(4, 2);
    for (int i = 0; i < 4; i++)
    {
        c->cd(i + 1);
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMyGraph(gx[i], "x Entry$", "eventID", 1, 4, 1);
        gx[i]->Draw();
        sprintf(buff, "Det%d:%d", i, idvecx[i].size());
        ll = DrawMyLatex(buff);
        ll->Draw();

        ll->Clear();
        c->cd(i + 5);
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMyGraph(gy[i], "y Entry$", "eventID", 1, 4, 1);
        sprintf(buff, "Det%d:%d", i, idvecy[i].size());
        ll = DrawMyLatex(buff);
        gy[i]->Draw();
        ll->Draw();
    }
    TString PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "IDvsEntry.png";
    c->SaveAs(PngName);
}

int CalculatePos(vector<AGETHit> (&buffer)[4])
{
    int counter = 0;
    for (int i = 0; i < 4; i++)
    {
        for (int i = 0; i < buffer[i].size(); i++)
        {
            counter++;
        }
    }
    return counter;
}
void ReadTrigger2Root(vector<TString> datList, TString fRawName, int force)
{
    if (datList.size() == 0)
    {
        cout << "Error!Dat size=0!" << endl;
        return;
    }
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fRawName << " is exist!" << endl;
            return;
        }
    }
    //TString fpath = GetPath(fRawName);
    //ofstream op(Form("%s/Trigger_event.dat", fpath.Data()));
    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
    {
        cout << "ERROR! " << fRawName << " cant open!" << endl;
        return;
    }
    //int ignore;
    int TrackID = 0;
    double ThTime = 0;
    double TOT = 0;
    TTree *tree = new TTree("tree", "");

    tree->Branch("TrackID", &TrackID, "TrackID/I");
    tree->Branch("LowThTime", &ThTime, "ThTime/D");
    //tree->Branch("LowTOT", &TOT, "TOT/D");

    counter = 0;
    int preTrackID = 0;
    int carry = 0;
    int breaknum = 0;
    int carrynum = 0;
    for (int i = 0; i < (int)datList.size(); i++)
    {

        FILE *infile;

        //sprintf(buff,"%s/CH%d_Data_%s_%s%s.txt",path,i+1,date,root,subtitle);
        infile = fopen(datList[i], "r");
        cout << "====>> Start to read File: " << datList[i] << endl;
        if (!infile)
        {
            cout << ">>> Invalid file!!! <<<" << endl;
            return;
        }

        while (!feof(infile))
        {
            TrackID = -999;
            ThTime = -999;
            fscanf(infile, "%d %lf", &TrackID, &ThTime);
            if (TrackID != -999)
            {

                //fscanf(infile[i], "%lf,%lf,%lf,%lf", &LowThTime[i], &HighThTime[i], &LowTOT[i], &HighTOT[i]);
                carry = IsCarry(TrackID, preTrackID);
                //op<<preTrackID<<"\t"<<TrackID;
                if (carry == 1e4)
                    breaknum++;
                if (carry == 65535)
                    carrynum++;
                preTrackID = TrackID;
                TrackID = breaknum * 1e4 + carrynum * 65535 + TrackID;
                //op<<"\t"<<TrackID<<endl;
#ifdef verbose
                cout << TrackID << "\t" << LowThTime << "\t" << HighThTime << "\t" << LowTOT << "\t" << HighTOT << endl;
                cout << TrackID << "\t" << LowThTime - HighThTime << endl;
#endif
                /*        
        if (preTrackID < 65536 && TrackID < preTrackID)
            mark++;
        preTrackID = TrackID;
        TrackID = mark * 65536 + TrackID;
#ifdef verbose
        cout << "#### before " << preTrackID << endl;
        cout << "#### after " << TrackID << endl;
#endif
*/
                tree->Fill();
                counter++;
            }
        }
        fclose(infile);
    }

    cout << "======>> The Event Number: " << counter << endl;
    cout << "====>> ReadFile progress is over , thank you ~" << endl;

    tree->Write();
    fFile->Close();
    //op.close();
}
void ReadFTOFData2Root(vector<TString> datList, TString fRawName, int force)
{
    if (datList.size() == 0)
    {
        cout << "Error!Dat size=0!" << endl;
        return;
    }
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fRawName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fRawName << " is exist!" << endl;
            return;
        }
    }

    TFile *fFile = new TFile(fRawName, "RECREATE");
    if (!fFile->IsOpen())
    {
        cout << "ERROR! " << fRawName << " cant open!" << endl;
        return;
    }

    //int fCH[] = {5, 6, 8, 12, 14};
    //int tureID[] = {2, 3, 100, 4, 1};
    //int ignore;
    int BoardID;
    int TrackID;
    int ChID;
    double LowThTime;
    double HighThTime;
    double LowTOT;
    double HighTOT;
    int Highflag; // if value=0, unvalid highthreshold information.
    TTree *tree = new TTree("tree", "");

    tree->Branch("TrackID", &TrackID, "TrackID/I");
    tree->Branch("BoardID", &BoardID, "BoardID/I");
    tree->Branch("ChID", &ChID, "ChID/I");
    tree->Branch("LowThTime", &LowThTime, "LowThTime/D");
    tree->Branch("HighThTime", &HighThTime, "HighThTime/D");
    tree->Branch("LowTOT", &LowTOT, "LowTOT/D");
    tree->Branch("HighTOT", &HighTOT, "HighTOT/D");
    tree->Branch("Highflag", &Highflag, "Highflag/I");

    counter = 0;
    int preTrackID = 0;
    int carry = 0;
    int breaknum = 0;
    int carrynum = 0;
    for (int i = 0; i < (int)datList.size(); i++)
    {
        FILE *infile;

        //sprintf(buff,"%s/CH%d_Data_%s_%s%s.txt",path,i+1,date,root,subtitle);
        infile = fopen(datList[i], "r");
        cout << "====>> Start to read File: " << datList[i] << endl;
        if (!infile)
        {
            cout << ">>> Invalid file!!! <<<" << endl;
            return;
        }

        while (!feof(infile))
        {
            ChID = -999;
            LowThTime = -999;
            HighThTime = -999;
            LowTOT = -999;
            HighTOT = -999;
            TrackID = -999;
            //fscanf(infile, "%d %d %lf %lf %lf %lf %d", &ignore, &ChID, &LowThTime, &HighThTime, &LowTOT, &HighTOT, &TrackID);
            fscanf(infile, "%d %d %lf %lf %lf %lf %d %d", &BoardID, &ChID, &LowThTime, &HighThTime, &LowTOT, &HighTOT, &TrackID, &Highflag);
            if (TrackID != -999)
            {
                carry = IsCarry(TrackID, preTrackID);
                //op<<preTrackID<<"\t"<<TrackID;
                if (carry == 1e4)
                    breaknum++;
                if (carry == 65535)
                    carrynum++;
                preTrackID = TrackID;
                TrackID = breaknum * 1e4 + carrynum * 65535 + TrackID;
                //fscanf(infile, "%d,%d,%lf,%lf,%lf,%lf,%d", &ignore, &ChID, &LowThTime, &HighThTime, &LowTOT, &HighTOT, &TrackID);
                //fscanf(infile[i], "%lf,%lf,%lf,%lf", &LowThTime[i], &HighThTime[i], &LowTOT[i], &HighTOT[i]);
#ifdef verbose
                cout << TrackID << "\t" << LowThTime << "\t" << HighThTime << "\t" << LowTOT << "\t" << HighTOT << endl;
                cout << TrackID << "\t" << LowThTime - HighThTime << endl;
#endif
                /*        
        if (preTrackID < 65536 && TrackID < preTrackID)
            mark++;
        preTrackID = TrackID;
        TrackID = mark * 65536 + TrackID;
#ifdef verbose
        cout << "#### before " << preTrackID << endl;
        cout << "#### after " << TrackID << endl;
#endif
*/
                tree->Fill();
                counter++;
            }
        }
        fclose(infile);
    }
    cout << "======>> The Event Number: " << counter << endl;
    cout << "====>> ReadFile progress is over , thank you ~" << endl;

    tree->Write();
    fFile->Close();
}
void CheckCrosstalk(vector<TimeData> &buffer)
{
    if(buffer.size()<1) return;
    vector<TimeData> Savedbuffer;
    double timemin = 10e3;
    double totmax = -999;
    for (int i = 0; i < buffer.size(); i++)
    {
        if (buffer[i].highthtime - buffer[i].lowthtime < 200 && buffer[i].lowtot>buffer[i].hightot)
            {

            Savedbuffer.push_back(buffer[i]);
            if(buffer[i].lowtot>totmax) totmax = buffer[i].lowtot;
            if(buffer[i].triggertimediff<timemin) timemin = buffer[i].triggertimediff;
            }
    }
if(Savedbuffer.size()<2) {

buffer.swap(Savedbuffer);
return;
}
// * time selection
    for (int j = 0; j < Savedbuffer.size() - 1; j++)
    {

        for (int i = 0; i < Savedbuffer.size() - j - 1; i++)
        {
            if (Savedbuffer[i].triggertimediff > Savedbuffer[i + 1].triggertimediff)
            {

                TimeData temp = Savedbuffer[i];
                Savedbuffer[i] = Savedbuffer[i + 1];
                Savedbuffer[i + 1] = temp;
            }
        }
    }
    for(int i=0; i<Savedbuffer.size(); i++)
    {
        if(Savedbuffer.back().triggertimediff-timemin>1e3) {
            vector<TimeData>::iterator itr = Savedbuffer.end()-1;
        Savedbuffer.erase(itr);
        }
        else break;
    }
if(Savedbuffer.size()<2) {

buffer.swap(Savedbuffer);
return;
}
// * tot selection
    for (int j = 0; j < Savedbuffer.size() - 1; j++)
    {

        for (int i = 0; i < Savedbuffer.size() - j - 1; i++)
        {
            if (Savedbuffer[i].lowtot > Savedbuffer[i + 1].lowtot)
            {

                TimeData temp = Savedbuffer[i];
                Savedbuffer[i] = Savedbuffer[i + 1];
                Savedbuffer[i + 1] = temp;
            }
        }
    }
    for(int i=0; i<Savedbuffer.size(); i++)
    {
        //cout<<"savedbuffer.size="<<Savedbuffer.size()<<endl;
        //cout<<"Savedbuffer.back().chid="<<Savedbuffer.back().chid<<endl;
        if(Savedbuffer.back().lowtot- FTOFtot[Savedbuffer.back().chid]< -2* FTOFtotsigma[Savedbuffer.back().chid]) {
            vector<TimeData>::iterator itr = Savedbuffer.end()-1;
        Savedbuffer.erase(itr);
        }
        else break;
    }
    //cout<<"program check"<<endl;
    if(Savedbuffer.size()<1) cout<<"error, please check your data"<<endl;
    else buffer.swap(Savedbuffer);
    //return;
}
void Checkrepetition(vector<TimeData> &buffer)
{
    vector<TimeData> Savedbuffer;
    vector<TimeData> event_per_ch[128];
    TimeData start;
    vector<double> trigger_time_diff;
    vector<double> low_tot;
    vector<int> chidvec;
    double bias = 9e5;
    bool flag = 0;
    int thesave = -1;
    int chkinds = 0; // the kinds of fired channel

    // only 1 event;
    if (buffer.size() < 2)
        return;

    for (int i = 0; i < buffer.size(); i++)
    {
        if (buffer.at(i).chid >= 128 || buffer.at(i).chid < 0)
            continue;
        event_per_ch[buffer.at(i).chid].push_back(buffer.at(i));
        trigger_time_diff.push_back(buffer.at(i).triggertimediff);
        low_tot.push_back(buffer.at(i).lowtot);
        chidvec.push_back(buffer.at(i).chid);
    }
    int min = TMath::MinElement(chidvec.size(), &chidvec.at(0));
    if (min < 0)
        min = 0;
    int max = TMath::MaxElement(chidvec.size(), &chidvec.at(0));
    if (max < 0)
        max = 0;
    for (int i = min; i < max + 1; i++)
    {
        if (event_per_ch[i].size() > 0)
        {
            chkinds++;
        }
    }
    // only 1 channel has events;
    // select the event that has the maximum tot
    if (chkinds < 2)
    {
        for (int i = min; i < max + 1; i++)
        {

            int maxtot = -999;
            thesave = -1;
            for (int j = 0; j < event_per_ch[i].size(); j++)
            {
                //cout << "Repeated ftof event is:" << endl;
                int temp = event_per_ch[i].at(j).lowtot;
                //cout << "    ch: " << event_per_ch[i].at(j).chid << ", lowtot: " << temp << endl;
                if (maxtot < temp && temp < 10e3)
                {
                    maxtot = temp;
                    thesave = j;
                }
            }
            if (thesave >= 0 && thesave < event_per_ch[i].size())
                Savedbuffer.push_back(event_per_ch[i].at(thesave));
        }
        buffer.swap(Savedbuffer);
        return;
    }

    // other cases
    // select the event that has the minimum bias;
    double mean_trigger_time_diff = TMath::Mean(trigger_time_diff.begin(), trigger_time_diff.end());
    for (int i = min; i < max + 1; i++)
    {

        thesave = -1;
        bias = 9e5;
        for (int j = 0; j < event_per_ch[i].size(); j++)
        {
            //cout << "Repeated ftof event is:" << endl;
            int temp = TMath::Abs(event_per_ch[i].at(j).triggertimediff - mean_trigger_time_diff);
            //cout << "    ch: " << event_per_ch[i].at(j).chid << ", bias: " << temp << endl;
            if (bias > temp)
            {
                bias = temp;
                thesave = j;
            }
        }
        if (thesave >= 0 && thesave < event_per_ch[i].size())
            Savedbuffer.push_back(event_per_ch[i].at(thesave));
    }
    buffer.swap(Savedbuffer);
    return;
}
int GetTriggerNum(TString fRawName2)
{
    int TriggerID;
    vector<int> triggeridvec;
    TFile *rawFile2 = new TFile(fRawName2, "read");
    TTree *tree2 = (TTree *)rawFile2->Get("tree");
    tree2->SetBranchAddress("TrackID", &TriggerID);
    int N = tree2->GetEntriesFast();

    int iEntry = 0;
    for (int i = 0, j = 0; i < N; i++)
    {
        DrawProcessbar(i, j, N);
        tree2->GetEntry(i);
        if (IsMatched2(TriggerID, triggeridvec, iEntry, 10))
            continue;
        triggeridvec.push_back(TriggerID);
    }
    RecordTriggers = triggeridvec.size();
    tree2->GetEntry(0);
    int start = TriggerID;
    tree2->GetEntry(N - 1);
    int end = TriggerID;
    NumofTriggers = end - start + 1;

    cout << " --> Provided Trigger Num: " << NumofTriggers << endl;
    cout << " --> Recorded Trigger Num: " << RecordTriggers << endl;
    rawFile2->Close();
    return RecordTriggers;
}
void AnalysisFTOF(TString fRawName1, TString fRawName2, TString fAnaName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fAnaName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fAnaName << " is exist!" << endl;
            return;
        }
    }

    ReadEleoffset();
    TFile *anaFile = new TFile(fAnaName, "recreate");
    if (!anaFile->IsOpen())
    {
        cout << "ERROR! " << fAnaName << " cant open!" << endl;
        return;
    }
    TimeData dethit;
    TTree *htree = new TTree("tree", "tree");
    htree->Branch("dethit", &dethit);
    int TrackID;
    int BoardID;
    int ChID;
    double LowThTime;
    double HighThTime;
    double LowTOT;
    double HighTOT;
    int Highflag;

    TFile *rawFile1 = new TFile(fRawName1, "read");
    TTree *tree1 = (TTree *)rawFile1->Get("tree");
    tree1->SetBranchAddress("TrackID", &TrackID);
    tree1->SetBranchAddress("BoardID", &BoardID);
    tree1->SetBranchAddress("ChID", &ChID);
    tree1->SetBranchAddress("LowThTime", &LowThTime);
    tree1->SetBranchAddress("HighThTime", &HighThTime);
    tree1->SetBranchAddress("LowTOT", &LowTOT);
    tree1->SetBranchAddress("HighTOT", &HighTOT);
    tree1->SetBranchAddress("Highflag", &Highflag);

    int TriggerID;
    double ThTime;
    double TOT;
    TFile *rawFile2 = new TFile(fRawName2, "read");
    TTree *tree2 = (TTree *)rawFile2->Get("tree");
    tree2->SetBranchAddress("TrackID", &TriggerID);
    tree2->SetBranchAddress("LowThTime", &ThTime);
    //tree2->SetBranchAddress("LowTOT", &TOT);
    vector<int> idvec;
    vector<int> identry;
    int N;
    N = tree1->GetEntriesFast();
    cout << " FTOF raw root entries: " << N << endl;
    for (int i = 0; i < N; i++)
    {
        tree1->GetEntry(i);
        if (LowTOT > 0)
        {
            idvec.push_back(TrackID);
            identry.push_back(i);
        }
    }

    N = tree2->GetEntriesFast();
    cout << " Trigger raw root entries: " << N << endl;
    NumofT0 = N;
    vector<int> triggeridvec;
    int iEntry = 0;
    int iEntrytrigger = 0;
    for (int i = 0, j = 0; i < N; i++)
    {

        DrawProcessbar(i, j, N);
        tree2->GetEntry(i);
        if (IsMatched2(TriggerID, triggeridvec, iEntrytrigger, 100))
            continue;
        triggeridvec.push_back(TriggerID);

        int nrange = 1000;
        int nentries = idvec.size();
        int start = (iEntry - nrange < 0) ? 0 : iEntry - nrange;
        int end = (iEntry + nrange > nentries) ? nentries : iEntry + nrange;
        for (int s = start; s < end; s++)
        {

            if (TriggerID == idvec[s])
            {
                iEntry = s;
                tree1->GetEntry(identry.at(s));
                dethit.event = TrackID;
                dethit.chid = GetDetChID(BoardID, ChID);
                dethit.lowtot = LowTOT;
                dethit.hightot = HighTOT;
                if (ChID >= 32)
                {

                    dethit.highthtime = HighThTime - T0EleTHoffset[BoardID] - T0Eleoffset[BoardID];
                    dethit.lowthtime = LowThTime - T0Eleoffset[BoardID];
                }
                else
                {

                    dethit.highthtime = HighThTime - FTOFEleTHoffset[BoardID * 32 + ChID] - FTOFEleoffset[BoardID * 32 + ChID];
                    dethit.lowthtime = LowThTime - FTOFEleoffset[BoardID * 32 + ChID];
                }
                dethit.triggertimediff = LowThTime - ThTime;
                //dethit.triggertot = TOT;

                htree->Fill();
            }
        }
    }

    anaFile->Flush();
    anaFile->WriteTObject(htree);
    anaFile->Close();
    rawFile1->Close();
    rawFile2->Close();
    cout << "--> Combined FTOF results has been stored to " << fAnaName << endl
         << endl;
}

void SortFTOF(TString fAnaName, TString fSorName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fSorName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fSorName << " is exist!" << endl;
            return;
        }
    }
    TFile *conFile = new TFile(fSorName, "recreate");
    if (!conFile->IsOpen())
    {
        cout << "ERROR! " << fSorName << " cant open!" << endl;
        return;
    }
    FTOFHit fhit;
    TTree *htree = new TTree("tree", "tree");
    htree->Branch("FTOFHit", &fhit);

    TimeData *dethit = new TimeData;

    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *ftoftree = (TTree *)anaFile->Get("tree");
    ftoftree->SetBranchAddress("dethit", &dethit);

    vector<int> ftofidvec;
    vector<int> ftofidvec_entrynum;
    vector<int> ftofidindexvec;
    vector<TimeData> ftofbuffer;
    int N;
    int idrange[2];
    N = ftoftree->GetEntriesFast();
    cout << " FTOF entries: " << N << endl;
    for (int i = 0; i < N; i++)
    {
        ftoftree->GetEntry(i);

        // pick up cosmicay signal
        if (dethit->triggertimediff > -700e3 && dethit->triggertimediff < -640e3)
        {
            ftofidvec.push_back(dethit->event);
            ftofidvec_entrynum.push_back(i);
        }
    }
    idrange[0] = TMath::MinElement(ftofidvec.size(), &ftofidvec.at(0));
    idrange[1] = TMath::MaxElement(ftofidvec.size(), &ftofidvec.at(0));
    for (int i = idrange[0]; i < idrange[1] + 1; i++)
    {
        ftofbuffer.clear();
        ftofidindexvec.clear();
        fhit = {0};
        if (SearchIDIndex(i, ftofidvec, ftofidvec_entrynum, ftofidindexvec))
        {

            for (int k = 0; k < ftofidindexvec.size(); k++)
            {
                int pos = ftofidindexvec.at(k);
                ftoftree->GetEntry(pos);
                ftofbuffer.push_back(*dethit);
            }
            Checkrepetition(ftofbuffer);
            fhit = Convert2T0DATA(ftofbuffer);
            htree->Fill();
        }
    }
    ValidT0 = htree->GetEntriesFast();
    conFile->Flush();
    conFile->WriteTObject(htree);
    conFile->Close();
    anaFile->Close();
    cout << "--> " << fAnaName << endl
         << " has been converted to:" << endl
         << "     -->" << fSorName << endl
         << endl;
}

void FilterFTOF(TString fAnaName, TString fConName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fConName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fConName << " is exist!" << endl;
            return;
        }
    }
    TFile *conFile = new TFile(fConName, "recreate");
    if (!conFile->IsOpen())
    {
        cout << "ERROR! " << fConName << " cant open!" << endl;
        return;
    }
    Readtot(100);

    FTOFHit fhit;
    TTree *htree = new TTree("tree", "tree");
    htree->Branch("FTOFhit", &fhit);

    TimeData *dethit = new TimeData;

    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *timetree = (TTree *)anaFile->Get("tree");
    timetree->SetBranchAddress("dethit", &dethit);

    // store event id
    vector<int> ftofidvec;
    vector<int> ftofidvec_entrynum;
    vector<int> ftofidindexvec;
    vector<TimeData> ftofbuffer;

    vector<int> T0idvec;
    vector<int> T0idvec_entrynum;
    vector<int> T0idindexvec;
    vector<TimeData> T0buffer;
    int N;
    int T0idrange[2];
    int FTOFidrange[2];
    N = timetree->GetEntriesFast();
    cout << " FTOF entries: " << N << endl;
    for (int i = 0; i < N; i++)
    {
        timetree->GetEntry(i);

        // pick up cosmicay signal
        if (dethit->chid >= 100)
        {
            if (dethit->triggertimediff < -380e3 && dethit->triggertimediff > -480e3)
            {
                T0idvec.push_back(dethit->event);
                T0idvec_entrynum.push_back(i);
            }
        }
        if (dethit->triggertimediff < -420e3 && dethit->triggertimediff > -500e3)
        //if (dethit->triggertimediff > -700e3 && dethit->triggertimediff < -640e3)
        {
            ftofidvec.push_back(dethit->event);
            ftofidvec_entrynum.push_back(i);
        }
    }
    if (T0idvec.size() >= 1)
    {

        T0idrange[0] = TMath::MinElement(T0idvec.size(), &T0idvec.at(0));
        T0idrange[1] = TMath::MaxElement(T0idvec.size(), &T0idvec.at(0));
    }
    else
    {
        cout << " T0 has no entries!" << endl;
    }
    if (ftofidvec.size() >= 1)
    {
        FTOFidrange[0] = TMath::MinElement(ftofidvec.size(), &ftofidvec.at(0));
        FTOFidrange[1] = TMath::MaxElement(ftofidvec.size(), &ftofidvec.at(0));
    }
    else
    {
        cout << " FTOF has no entries!" << endl;
    }
    int min = T0idrange[0] < FTOFidrange[0] ? T0idrange[0] : FTOFidrange[0];
    int max = T0idrange[1] < FTOFidrange[1] ? T0idrange[1] : FTOFidrange[1];
    int FTOFentry = 0;
    int T0entry = 0;
    for (int i = min, j = 0; i < max + 1; i++)
    //for (int i = 7880; i < 15000; i++)
    {
        DrawProcessbar(i - min, j, max + 1 - min);
        ftofbuffer.clear();
        ftofidindexvec.clear();
        T0buffer.clear();
        T0idindexvec.clear();
        fhit.initial();
        //if (SearchIDIndex(i, ftofidvec, ftofidvec_entrynum, ftofidindexvec))
        if (SearchIDIndex2(i, ftofidvec, ftofidvec_entrynum, ftofidindexvec, FTOFentry, 500))
        {

            for (int k = 0; k < ftofidindexvec.size(); k++)
            {
                int pos = ftofidindexvec.at(k);
                timetree->GetEntry(pos);
                ftofbuffer.push_back(*dethit);
            }
            Checkrepetition(ftofbuffer);
        }
        //cout<<"ftofbuffer.size="<<ftofbuffer.size()<<endl;
        CheckCrosstalk(ftofbuffer);
        //cout<<"check crosstalk, N="<<i<<endl;
        //if (SearchIDIndex(i, T0idvec, T0idvec_entrynum, T0idindexvec))
        if (SearchIDIndex2(i, T0idvec, T0idvec_entrynum, T0idindexvec, T0entry, 500))
        {

            for (int k = 0; k < T0idindexvec.size(); k++)
            {
                int pos = T0idindexvec.at(k);
                timetree->GetEntry(pos);
                T0buffer.push_back(*dethit);
            }
            Checkrepetition(T0buffer);
        }
        CheckCrosstalk(T0buffer);
        if (ftofbuffer.size() > 0 || T0buffer.size() > 0)
        {

            T0buffer.insert(T0buffer.end(), ftofbuffer.begin(), ftofbuffer.end());
            //cout<<"after filter, fires= "<<T0buffer.size()<<endl;
            fhit = Convert2T0DATA(T0buffer);
            //cout<<"> progress check < "<<i<<endl;
            htree->Fill();
        }
    }
    ValidT0 = htree->GetEntriesFast();
    conFile->Flush();
    conFile->WriteTObject(htree);
    conFile->Close();
    anaFile->Close();
    cout << "--> " << fAnaName << endl
         << " has been converted to:" << endl
         << "     -->" << fConName << endl
         << endl;
}

void Match2(TString fName1, TString fName2, TString fMatName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fMatName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fMatName << " is exist!" << endl;
            return;
        }
    }

    bool Trackerflag;
    bool T0flag;
    //* FTOF data
    FTOFHit *FTOFhit = new FTOFHit;
    TFile *t0file;
    TTree *t0tree;
    if (fName1.CompareTo(""))
    {
        t0file = new TFile(fName1, "read");
        t0tree = (TTree *)t0file->Get("tree");
        t0tree->SetBranchAddress("FTOFhit", &FTOFhit);
        T0flag = 1;
    }

    //* Tracker data
    MMHit *MMhit = new MMHit;
    TFile *trackerfile;
    TTree *trackertree[NumMM];
    if (fName2.CompareTo(""))
    {
        trackerfile = new TFile(fName2, "read");
        Trackerflag = 1;
    }
    if (!T0flag && !Trackerflag)
        return;
    int allmin = 9999999;
    int allmax = 0;
    for (int s = 0; s < NumMM; s++)
    {
        trackertree[s] = (TTree *)trackerfile->Get(Form("tree%d", s));
        trackertree[s]->SetBranchAddress("MMhit", &MMhit);
        cout << "Entry Detector " << s << " ..." << endl;
        int nEntries = trackertree[s]->GetEntriesFast();
        trackertree[s]->GetEntry(0);
        int min = MMhit->event;
        allmin = (min < allmin) ? min : allmin;
        trackertree[s]->GetEntry(nEntries - 1);
        int max = MMhit->event;
        allmax = (max > allmax) ? max : allmax;
        cout << "--> "
             << "Tracker" << s << ", Events: " << nEntries << ", tigger range: " << min << " " << max << endl;
    }

    cout << "--> Searching trigger from " << allmin << " to " << allmax << endl;

    // create new tree;
    TFile *comFile = new TFile(fMatName, "recreate");
    if (!comFile->IsOpen())
    {
        cout << "ERROR! " << fMatName << " cant open!" << endl;
        return;
    }
    TTree *tree = new TTree("tree", "tree");
    CRHit fHit;
    tree->Branch("CRhit", &fHit);
    int iEntry[NumMM] = {0};
    for (int trigid = allmin; trigid < allmax && trigid >= 0 && trigid < 100e4; trigid++)
    {
        if (trigid % 100 == 0)
            cout << "--> Searching trigger id = " << trigid << endl;
        fHit.initial();
        fHit.event = trigid;
        //for tracker
        for (int s = 0; s < NumMM; s++)
        {
            int nrange = 100;
            int nentries = trackertree[s]->GetEntriesFast();
            int start = (iEntry[s] - nrange < 0) ? 0 : iEntry[s];
            int end = (iEntry[s] + nrange > nentries) ? nentries : iEntry[s] + nrange;
            MMHit MMtemp;
            vector<MMHitData> MMdatatemp[2];
            for (int i = start; i < end; i++)
            {
                MMhit->initial();
                trackertree[s]->GetEntry(i);
                if (MMhit->event == trigid)
                {

                    iEntry[s] = i;
                    //PushBack(fEventList->hit, fEvent);
                    if (MMhit->xyflag == 1)
                    {

                        cout << "MM" << s << ",x have cluster" << endl;
                        MMtemp = *MMhit;
                        MMdatatemp[0] = MMhit->MMchn[0];
                    }

                    if (MMhit->xyflag == 2)
                    {
                        cout << "MM" << s << ",y have cluster" << endl;
                        MMtemp = *MMhit;
                        MMdatatemp[1] = MMhit->MMchn[1];
                    }
                    //break;
                }
            }
            if (MMdatatemp[1].size() && MMdatatemp[0].size())
            {

                MMtemp.xyflag = 3;
                MMtemp.MMchn[0] = MMdatatemp[0];
                MMtemp.MMchn[1] = MMdatatemp[1];
            }
            if (MMdatatemp[1].size() || MMdatatemp[0].size())
            {

                fHit.Tracker.push_back(MMtemp);
                /*
                cout << "xyflag:" << MMtemp.xyflag;
                for (int ii = 0; ii < MMtemp.MMchn[0].size(); ii++)
                {
                    cout << " x chn:" << MMtemp.MMchn[0][ii].clustermean[0];
                }
                for (int ii = 0; ii < MMtemp.MMchn[1].size(); ii++)
                {
                    cout << " y chn:" << MMtemp.MMchn[1][ii].clustermean[1];
                }
                cout << endl;
                */
            }
        }

        // T0?
        if (T0flag)
        {
            //
        }
        // FTOF?

        if (fHit.Tracker.size())
            tree->Fill();
    }
    cout << "--> " << tree->GetEntries() << " events are stored." << endl;
    comFile->WriteTObject(tree);
    comFile->Close();
    //if(t0file->IsOpen()) t0file->Close();
    //trackerfile->Close();
}

void AnalysisTrack(vector<int> uselist, int target, vector<MMHit> MMs, int xyflag)
{
    TGraphErrors *g1 = new TGraphErrors(); //方法1: 只要是好的cluster，就使用这个信息
    TGraphErrors *g2 = new TGraphErrors(); //方法2: 只用q最大的那个好的cluster
    int ipoint1 = 0;
    int ipoint2 = 0;
    for (int i = 0; i < (int)uselist.size(); i++)
    {
        int id = -1;
        int detID = uselist[i];
        for (int j = 0; j < MMs.size(); j++)
            if (MMs[j].MMID == detID)
                id = j;
        if (id == -1)
            continue;

        double qmax = -1;
        double mmax = -1;
        vector<MMHitData> cluster = MMs[id].MMchn[xyflag];
        vector<MMHitData> validcluster;

        double Nbranch = cluster.size();
        for (int jj = 0; jj < cluster.size(); jj++)
        {
            if (!IsTrackClusterEffective(cluster.at(jj)))
                continue;
            validcluster.push_back(cluster.at(jj));
        }
        // number of cluster <=3
        if (validcluster.size() > 3)
            continue;
        for (int jj = 0; jj < validcluster.size(); jj++)
        {
            double xy = CalTrackPosition(id, xyflag, validcluster.at(jj).clustermean[xyflag]);
            double z = CalTrackPosition(id, Z);
            g1->SetPoint(ipoint1++, z, xy);

            mmax = (qmax < validcluster.at(jj).q) ? validcluster.at(jj).clustermean[xyflag] : mmax;
            qmax = (qmax < validcluster.at(jj).q) ? validcluster.at(jj).q : qmax;
        }
        if (mmax != -1)
        {
            double xy = CalTrackPosition(id, xyflag, mmax);
            double z = CalTrackPosition(id, Z);
            g2->SetPoint(ipoint2++, z, xy);
        }
    }
    //拟合
    if (g1->GetN() < 2)
        return;
    if (g2->GetN() < 2)
        return;
    g1->Fit("pol1", "q");
    g2->Fit("pol1", "q");
    //target上击中的期望值
    double exp1 = g1->GetFunction("pol1")->Eval(CalTrackPosition(target, Z));
    double exp2 = g2->GetFunction("pol1")->Eval(CalTrackPosition(target, Z));

    //找出target在detlist里的编号
    int id = -1;
    for (int j = 0; j < MMs.size(); j++)
        if (MMs[j].MMID == target)
            id = j;

    //如果target没有击中，则填默认值
    if (id == -1)
    {
        fDeviation[target][xyflag]->Fill(-999);
        return;
    }
    vector<MMHitData> cluster = MMs[id].MMchn[xyflag];
    //循环target的所有击中，求残差
    for (int j = 0; j < (int)cluster.size(); j++)
    {
        if (!IsTrackClusterEffective(cluster.at(j)))
            continue;

        double xy = CalTrackPosition(target, xyflag, cluster.at(j).clustermean[xyflag]);
        fDeviation[target][xyflag]->Fill(xy - exp2);
    }

    delete g1;
    delete g2;
}

void FillTrackHistogram(vector<MMHit> detlist, int detID)
{

    for (int i = 0; i < (int)detlist.size(); i++)
    {
        if (detlist[i].MMID != detID)
            continue;
        MMHit detector = detlist[i];
        // number of cluster 分布
        for (int xy = 0; xy < 2; xy++)
        {
            fNCluster[detID][xy]->Fill(detector.MMchn[xy].size());
            for (int j = 0; j < (int)detector.MMchn[xy].size(); j++)
            {

                fCluster[detID][xy]->Fill(detector.MMchn[xy][j].stripnum);

                fTotCharge[detID][xy]->Fill(detector.MMchn[xy][j].q);
            }
        }
        fNCluster2D[detID]->Fill(detector.MMchn[X].size(), detector.MMchn[Y].size());

        //fhitmap分布
        double xmean, xrms, xq;
        double ymean, yrms, yq;
        double xqmax = -1;
        double xmmax = -1;
        double yqmax = -1;
        double ymmax = -1;
        for (int ii = 0; ii < (int)detector.MMchn[X].size(); ii++)
            for (int jj = 0; jj < (int)detector.MMchn[Y].size(); jj++)
            {

                fhitmap[detID][0]->Fill(detector.MMchn[X][ii].clustermean[0], detector.MMchn[Y][jj].clustermean[1]);
                if (!IsTrackClusterEffective(detector.MMchn[X][ii]))
                    continue;
                if (!IsTrackClusterEffective(detector.MMchn[Y][jj]))
                    continue;

                fhitmap[detID][1]->Fill(detector.MMchn[X][ii].clustermean[0], detector.MMchn[Y][jj].clustermean[1]);

                //在束流实验室系下的坐标在这里写出来
                //cout << " x strip: " << detector.MMchn[X][ii].clustermean[0] << ", y strip: " << detector.MMchn[Y][jj].clustermean[1] << endl;
                double x = CalTrackPosition(detID, X, detector.MMchn[X][ii].clustermean[0]);
                double y = CalTrackPosition(detID, Y, detector.MMchn[Y][jj].clustermean[1]);
                //cout << " x pos: " << x << ", y pos: " << y << endl;
                fhitmapReal[detID][0]->Fill(x, y);
                xmmax = (xqmax < detector.MMchn[X][ii].q) ? detector.MMchn[X][ii].clustermean[0] : xmmax;
                xqmax = (xqmax < detector.MMchn[X][ii].q) ? detector.MMchn[X][ii].q : xqmax;
                ymmax = (yqmax < detector.MMchn[Y][jj].q) ? detector.MMchn[Y][jj].clustermean[1] : ymmax;
                yqmax = (yqmax < detector.MMchn[Y][jj].q) ? detector.MMchn[Y][jj].q : yqmax;
            }
        if (xmmax != -1 || ymmax != -1)
        {
            //cout << "way2: x strip: " << xmmax << ", y strip: " << ymmax << endl;
            double x = CalTrackPosition(detID, X, xmmax);
            double y = CalTrackPosition(detID, Y, ymmax);
            //cout << "way2: x pos: " << x << ", y pos: " << y << endl;
            fhitmapReal[detID][1]->Fill(x, y);
        }
        //cout << detName[detlist[i].MMID] << " hist has been filled " << endl;
    }
}

bool IsTrackClusterEffective(MMHitData cluster)
{

    //判断是否是放电信号或噪声信号

    if (cluster.stripnum < 2)
        //if (cluster.q < 5 || cluster.stripnum < 2 || cluster.q > 500 || cluster.stripnum > 10)
        return false;
    return true;
}
double CalTrackPosition(int detID, int flag, double channel)
{
    double position = 0;
    if (detID == FTOFid)
        return position;
    if (detID == T0id)
        return position;

    if (detID == MM0id || detID == MM1id || detID == MM2id || detID == MM3id)
    {
        if (flag == X)
            position = Xoff[detID] + (384 / 2. - channel + 1 / 2.) * 0.4;
        if (flag == Y)
            position = Yoff[detID] + (channel - 384 / 2. + 1 / 2.) * 0.4;
        if (flag == Z)
            position = Zpos[detID];
    }
    return position;
}
void Definehist()
{
    //--------------------------
    // 定义histogram & graph
    int nbin;

    /*
    ftmp1 = new TH1F("ftmp1", "tmp", 128, 0, 128);
    ftmp2 = new TH2F("ftmp2", "tmp", 32, 0, 32, 32, 0, 32);
    //光子重建相关
    fRec = new TH1F("fRec", "reconstructed cherenkov angle by single cluster", 100, 0.5, 1.5);           //3.1415926/2);
    fRec2 = new TH1F("fRec2", "reconstructed mean cherenkov angle by multiple clusters", 100, 0.5, 1.5); //3.1415926/2);
    fNRec = new TH1F("fNRec", "number reconstructed cherenkov photons", 10, 0, 10);                      //3.1415926/2);
    fQphoton = new TH1F("fQphoton", "charge for photon signal", 256, 0, 4096 * 5);
    fTphoton = new TH1F("fTphoton", "time for photon signal", 100, 200, 300);
    fQTphoton = new TH2F("fQTphoton", "Q vs T for photons", 256, 0, 4096, 100, 200, 300);
    fTotQphoton = new TH1F("fTotQphoton", "total photon charges", 256, 0, 4096);
    fNClusterPhoton = new TH1F("fNClusterPhoton", "number of cluster for photons", 10, 0, 10);
    fClusterPhoton = new TH1F("fClusterPhoton", "cluster size for photons", 10, 0, 10);

    gRec = new TGraph2D();
    gRec->GetXaxis()->SetTitle("X");
    gRec->GetYaxis()->SetTitle("Y");
    fFootmap = new TGraph(); //new TH2F("fFootmap", "foot map for RICH",  40, -40 / 2 * 5, 40 / 2 * 5, 40, -40 / 2 * 5, 40 / 2 * 5);
*/
    fHitmap = new TH2F("fHitmap", "Real Position hit map", 40, -40 / 2 * 5, 40 / 2 * 5, 40, -40 / 2 * 5, 40 / 2 * 5);
    fHitmap->GetXaxis()->SetTitle("X");
    fHitmap->GetYaxis()->SetTitle("Y");
    //事例填图
    for (int i = 0; i < NumMM; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            /*
            */
            nbin = 384;
            fhitmap[i][j] = new TH2F(Form("fhitmap%d_%d", i, j), Form("hit map for %s", detName[i]), nbin, 0, nbin, nbin, 0, nbin);
            fhitmap[i][j]->GetXaxis()->SetTitle("X");
            fhitmap[i][j]->GetYaxis()->SetTitle("Y");

            //nbin = (i == T06id) ? 25 : 5;
            nbin = 20;
            fDeviation[i][j] = new TH1F(Form("fdev%d_%d", i, j), Form("Deviation of Track %s for %s", detName[i], xyName[j]), 100, -1 * nbin, nbin);
            fDeviation[i][j]->SetXTitle(Form("%s(mm)", xyName[j]));
            fDeviation[i][j]->SetYTitle("Entries");

            fNCluster[i][j] = new TH1F(Form("fnclu%d_%d", i, j), Form("Number of cluster of %s for %s", detName[i], xyName[j]), 10, 0, 10);
            fCluster[i][j] = new TH1F(Form("fclu%d_%d", i, j), Form("Cluster size of %s for %s", detName[i], xyName[j]), 10, 0, 10);
            fhitmapReal[i][j] = new TH2F(Form("fhitmapReal%d_%d", i, j), Form("Real Position hit map for %s", detName[i]), 800, -800 / 2 * 0.4, 800 / 2 * 0.4, 800, -800 / 2 * 0.4, 800 / 2 * 0.4);

            fTotCharge[i][j] = new TH1F(Form("ftotQ%d_%d", i, j), Form("Total charge of %s for %s", detName[i], xyName[j]), 1e3, 0, 1000);
            fTotCharge[i][j]->GetXaxis()->SetTitle("Q");
            fTotCharge[i][j]->GetYaxis()->SetTitle("Entries");
        }
        //nbin = (i == 0) ? 32 : 128;
        fNCluster2D[i] = new TH2F(Form("fnclu2%d", i), Form("Number of cluster of %s", detName[i]), 10, 0, 10, 10, 0, 10);
        /*
        fPeakingTime[i] = new TH1F(Form("fT%d", i), Form("Peaking time of %s", detName[i]), 100, 0, 3000);
        fPeakingTime[i]->GetXaxis()->SetTitle("T");
        fPeakingTime[i]->GetYaxis()->SetTitle("Entries");



        fQT[i] = new TH2F(Form("fQT%d", i), Form("Q vs. T of %s", detName[i]), 256, 0, 4096, 100, 200, 300);
        fQT[i]->GetXaxis()->SetTitle("Q");
        fQT[i]->GetYaxis()->SetTitle("T");
*/
    }
}
void DrawTracker(TString path, int id)
{
    TCanvas *cc = new TCanvas(Form("c%d", id));
    cc->SetTitle(detName[id]);
    cc->Clear();
    cc->Divide(4, 4);
    cc->cd(1);
    fhitmap[id][0]->Draw("colz");
    cc->cd(2);
    fNCluster[id][0]->Draw("");
    cc->cd(3);
    fCluster[id][0]->Draw("");
    cc->cd(4);
    fTotCharge[id][0]->Draw("");
    cc->cd(5);
    fhitmap[id][1]->Draw("colz");
    cc->cd(6);
    fNCluster[id][1]->Draw("");
    cc->cd(7);
    fCluster[id][1]->Draw("");
    cc->cd(8);
    fTotCharge[id][1]->Draw("");
    cc->cd(9);
    fhitmapReal[id][0]->Draw("colz");
    cc->cd(10);
    fhitmapReal[id][0]->ProjectionX()->Draw();
    //fPeakingTime[id]->Draw("");
    cc->cd(11);
    fhitmapReal[id][0]->ProjectionY()->Draw();
    //fCharge[id]->Draw("");
    cc->cd(12);
    fNCluster2D[id]->Draw("boxtext");
    //fQT[id]->Draw("colz");
    cc->cd(13);
    fhitmapReal[id][1]->Draw("colz");
    cc->cd(14);
    fhitmapReal[id][1]->ProjectionX()->Draw();
    cc->cd(15);
    fhitmapReal[id][1]->ProjectionY()->Draw();
    cc->SaveAs(Form("%s/%strack.png", path.Data(), detName[id]));
}

void DrawDeviation(TString cname, vector<int> uselist)
{
    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle("Deviation distribution for " + fName);
    cc->Clear();
    cc->Divide(uselist.size(), 2);
    for (int i = 0; i < (int)uselist.size(); i++)
    {
        cc->cd(1 + i);
        if (fDeviation[uselist[i]][0]->Integral())
        {
            int fitRL = fDeviation[uselist[i]][0]->GetMean() - fDeviation[uselist[i]][0]->GetRMS();
            int fitRR = fDeviation[uselist[i]][0]->GetMean() + fDeviation[uselist[i]][0]->GetRMS();
            fDeviation[uselist[i]][0]->Fit("gaus", "q", "", fitRL, fitRR);
        }
        cc->cd(1 + i + uselist.size());
        if (fDeviation[uselist[i]][1]->Integral())
        {
            int fitRL = fDeviation[uselist[i]][1]->GetMean() - fDeviation[uselist[i]][1]->GetRMS();
            int fitRR = fDeviation[uselist[i]][1]->GetMean() + fDeviation[uselist[i]][1]->GetRMS();
            fDeviation[uselist[i]][1]->Fit("gaus", "q", "", fitRL, fitRL);
        }

        cout << detName[uselist[i]]
             << " efficiency_x = " << fDeviation[uselist[i]][0]->Integral(1, 100) / fDeviation[uselist[i]][0]->GetEntries()
             << " efficiency_y = " << fDeviation[uselist[i]][1]->Integral(1, 100) / fDeviation[uselist[i]][1]->GetEntries() << endl;
    }
}

//---------------
// 读取offset
void ReadOffset()
{
    for (int i = 0; i < (int)offset.size(); i++)
    {
        if (fName == TString(Form("RUN%02d", (int)offset[i][0])))
        {
            cout << "==> Reading offset from : " << offset[i][0] << ". " << endl;
            for (int j = 0; j < NumMM; j++)
            {
                Xoff[j] = offset[i + 0][j + 1];
                Yoff[j] = offset[i + 1][j + 1];
                Zpos[j] = offset[i + 2][j + 1];
            }
            angle = offset[i + 3][1];
            return;
        }
    }

    cout << "Warning: can't find the offset for " << fName << ". Please check." << endl;
}
//---------------
// 输出信息
void PrintOffset(int xyz)
{
    if (xyz == ANGLE)
    {
        cout << "{" << fName.ReplaceAll("RUN", "") << ", " << angle << "}, // angle" << endl;
        return;
    }

    cout.precision(4);
    cout << "{" << fName.ReplaceAll("RUN", "") << ", ";
    for (int i = 0; i < NumMM; i++) //FG125111
    {
        if (xyz == X)
        {

            //if (i != T02id && i != T03id && fDeviation[i][0]->GetFunction("gaus") != 0) //固定前两个tracker去标定后面几个tracker的位置
            if (fDeviation[i][0]->Integral())
                cout << -1 * fDeviation[i][0]->GetFunction("gaus")->GetParameter(1) + Xoff[i];
            //else
        }

        if (xyz == Y)
        {

            //if (i != T02id && i != T03id && fDeviation[i][1]->GetFunction("gaus") != 0)
            if (fDeviation[i][1]->Integral())
                cout << -1 * fDeviation[i][1]->GetFunction("gaus")->GetParameter(1) + Yoff[i];
            //else
        }
        if (xyz == Z)
            cout << Zpos[i];
        if (i != NumMM - 1)
            cout << ", ";
    }
    cout << "}, // " << xyName[xyz] << endl;
    cout.precision(6);
}
void PrintOffset2(int xyz)
{
    if (xyz == ANGLE)
    {
        cout << "{" << fName.ReplaceAll("RUN", "") << ", " << angle << "}, // angle" << endl;
        return;
    }

    cout.precision(4);
    cout << "{" << fName.ReplaceAll("RUN", "") << ", ";
    for (int i = 0; i < NumMM; i++) //FG125111
    {
        if (xyz == X)
        {

            //if (i != T02id && i != T03id && fDeviation[i][0]->GetFunction("gaus") != 0) //固定前两个tracker去标定后面几个tracker的位置
            if (fhitmapReal[i][0]->Integral())
                cout << -1 * fhitmapReal[i][0]->GetMean(1) + Xoff[i];
        }

        if (xyz == Y)
        {

            //if (i != T02id && i != T03id && fDeviation[i][1]->GetFunction("gaus") != 0)
            if (fhitmapReal[i][1]->Integral())
                cout << -1 * fhitmapReal[i][1]->GetMean(1) + Yoff[i];
        }
        if (xyz == Z)
            cout << Zpos[i];
        if (i != NumMM - 1)
            cout << ", ";
    }
    cout << "}, // " << xyName[xyz] << endl;
    cout.precision(6);
}
void PrintInfo()
{
    cout << "\n\n----------------Offset for " << fName << "-------------------\n"
         << endl;
    cout << "use hist fdeviation:";
    PrintOffset(X);
    PrintOffset(Y);
    PrintOffset(Z);
    PrintOffset(ANGLE);

    cout << "use hist fhitmapReal:";
    PrintOffset2(X);
    PrintOffset2(Y);
    PrintOffset2(Z);
    PrintOffset2(ANGLE);
}

void Matchall(TString fName1, TString fName2, TString fMatName, int force)
{
    TString fpath = GetPath(fMatName);
    ofstream idT0(Form("%s/idT0.dat", fpath.Data()));
    ofstream idTrackerx;
    ofstream idTrackery;
    ofstream idMatched(Form("%s/idMatched.dat", fpath.Data()));

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fMatName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fMatName << " is exist!" << endl;
            return;
        }
    }

    TFile *comFile = new TFile(fMatName, "recreate");
    if (!comFile->IsOpen())
    {
        cout << "ERROR! " << fMatName << " cant open!" << endl;
        return;
    }

    FTOFHit fhit;
    FTOFHit *FTOFhit = new FTOFHit;
    AGETHit *trackerhit = new AGETHit;
    // create new tree;
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("FTOFHit", &fhit);
    //* T0 data
    TFile *t0file = new TFile(fName1, "read");
    TTree *t0tree = (TTree *)t0file->Get("tree");
    t0tree->SetBranchAddress("FTOFhit", &FTOFhit);

    //* Tracker data
    TFile *trackerfile = new TFile(fName2, "read");
    TTree *trackertree[4];
    vector<AGETHit> trackerbuffer[4];
    vector<int> idvecx[4];
    vector<int> idvecx_entrynum[4];
    vector<int> idvecy[4];
    vector<int> idvecy_entrynum[4];
    int trackercounter = 0;
    for (int s = 0; s < 4; s++)
    {
        idTrackerx.open(Form("%s/idTrackerx_%d.dat", fpath.Data(), s));
        idTrackery.open(Form("%s/idTrackery_%d.dat", fpath.Data(), s));

        trackertree[s] = (TTree *)trackerfile->Get(Form("tree%d", s));
        trackertree[s]->SetBranchAddress("dethit", &trackerhit);
        cout << "Entry Detector " << s << " ..." << endl;

        idvecx[s].clear();
        idvecy[s].clear();

        int nentries = trackertree[s]->GetEntriesFast();
        cout << "Entries is " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {
            Long64_t ii = trackertree[s]->LoadTree(i);
            if (ii < 0)
                break;

            trackertree[s]->GetEntry(ii);
            //cout<<"event id: "<<trackerhit->event<<endl;
            // x axis has a hit
            if (trackerhit->wave_max > 150 && trackerhit->q > 6)
            //if (trackerhit->wave_max > 7 * trackerhit->ped_rms && trackerhit->q > 5)
            //if (1)
            {

                if (trackerhit->chn[0] != -999)
                {
                    if (idvecx[s].size() < 1)
                    {

                        idvecx[s].push_back(trackerhit->event);
                        idTrackerx << idvecx[s].back() << endl;
                    }
                    if (idvecx[s].back() != trackerhit->event)
                    {

                        idvecx[s].push_back(trackerhit->event);
                        idTrackerx << idvecx[s].back() << endl;
                    }
                }

                // y axis has a hit
                if (trackerhit->chn[1] != -999)
                {
                    if (idvecy[s].size() < 1)
                    {

                        idvecy[s].push_back(trackerhit->event);
                        idTrackery << idvecy[s].back() << endl;
                    }
                    if (idvecy[s].back() != trackerhit->event)
                    {

                        idvecy[s].push_back(trackerhit->event);
                        idTrackery << idvecy[s].back() << endl;
                    }
                }
            }
        }
        idTrackerx.close();
        idTrackery.close();
    }
    //sort(idvecx[0].begin(), idvecx[0].end());
    //sort(idvecy[0].begin(), idvecy[0].end());

    // Tracker X Y match

    for (int s = 0; s < 4; s++)
    {
        int trackerxwaitmatched = idvecx[s].size();
        int trackerywaitmatched = idvecy[s].size();
        int matchedtracker = 0;
        for (int i = 0; i < trackerxwaitmatched; i++)
        {
            //oplog << "Search event id: " << FTOFHit->event << endl;
            //cout << "Search event id: " << FTOFHit->event << endl;
            //if (IsMatched(FTOFHit->event, idvecx[s]))
            //    matchedcounter++;

            if (IsMatched(idvecx[s].at(i), idvecy[s]))
                matchedtracker++;

            //if (1)
            // ##TODO
        }
        cout << "--> Det" << s << ", y hit: " << trackerywaitmatched << ", x hit: " << trackerxwaitmatched << endl
             << "    |-->the matched: " << matchedtracker << endl
             << endl;
    }

    int waitmatched = t0tree->GetEntriesFast();
    int matchedx = 0;
    int matchedy = 0;
    int matchedtotal = 0;
    int matchedor = 0;
    for (int i = 0; i < waitmatched; i++)
    {
        t0tree->GetEntry(i);
        //oplog << "Search event id: " << FTOFHit->event << endl;
        //cout << "Search event id: " << FTOFHit->event << endl;
        int matchedcounterx = 0;
        int matchedcountery = 0;
        int matchedcountertotal = 0;
        int xflag = 0;
        int yflag = 0;
        idT0 << FTOFhit->event << endl;
        for (int s = 0; s < 4; s++)
        {
            xflag = 0;
            yflag = 0;
            if (IsMatched(FTOFhit->event, idvecx[s]))
            {
                matchedcounterx++;
                xflag = 1;
            }
            if (IsMatched(FTOFhit->event, idvecy[s]))
            {
                matchedcountery++;
                yflag = 1;
            }
            if (xflag && yflag)
                matchedcountertotal++;
        }
        //if (1)
        // ##TODO
        if (matchedcounterx > 0)
            matchedx++;
        if (matchedcountery > 0)
            matchedy++;
        if (matchedcountery + matchedcounterx > 0)
            matchedor++;
        if (matchedcountertotal > 0)
        {
            matchedtotal++;
            fhit = *FTOFhit;
            tree->Fill();
            idMatched << FTOFhit->event << endl;
        }
    }
    comFile->Flush();
    comFile->WriteTObject(tree);
    comFile->Close();
    t0file->Close();
    trackerfile->Close();
    idT0.close();
    idMatched.close();
    cout << "-->  Matched data has been saved :" << endl
         << fMatName << endl
         << "    |-->  The T0 wait-matched Numbers : " << waitmatched << endl
         << "    |-->  The matched x Numbers : " << matchedx << ", Eff : " << (double)matchedx / waitmatched << endl
         << "    |-->  The matched y Numbers : " << matchedy << ", Eff : " << (double)matchedy / waitmatched << endl
         << "    |-->  The matched x&y Numbers : " << matchedtotal << ", Eff : " << (double)matchedtotal / waitmatched << endl
         << "    |-->  The matched x|y Numbers : " << matchedor << ", Eff : " << (double)matchedor / waitmatched << endl
         << endl;
}

void TrackerEff(TString fAnaName)
{
    FileStat_t fStat;
    gSystem->GetPathInfo(fAnaName, fStat);
    if (fStat.fSize == 0)
    {
        cout << "ERROR! " << fAnaName << " isn't exist!" << endl;
        return;
    }
    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[4];
    AGETHit *dethit;
    dethit = new AGETHit;
    TString fpath = GetPath(fAnaName);
    ofstream op(Form("%s/Tracker_effinfo.dat", fpath.Data()));

    TH1I *hid = new TH1I("hid", "hid", 100e4, 1, 100e4 + 1);
    TH1I *hXid[4];
    TH1I *hYid[4];
    TH1I *hXHit[4];
    TH1I *hYHit[4];
    int chNo = 64; // 64*6=384
    vector<int> xid[4];
    vector<int> yid[4];
    hid->Clear();
    for (int s = 0; s < 4; s++)
    {
        hXid[s] = new TH1I(Form("hXid%d", s), Form("hXid%d", s), 100e4, 1, 100e4 + 1);
        hYid[s] = new TH1I(Form("hYid%d", s), Form("hYid%d", s), 100e4, 1, 100e4 + 1);
        hXHit[s] = new TH1I(Form("HitX%d", s), "X Strip hits for single event", chNo, 0, chNo);
        hYHit[s] = new TH1I(Form("HitY%d", s), "Y Strip hits for single event", chNo, 0, chNo);
        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        tree[s]->SetBranchAddress("dethit", &dethit);
        int nentries = tree[s]->GetEntriesFast();
        //nentries =1e4;
        cout << " --> nentries= " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {

            Long64_t ii = tree[s]->LoadTree(i);
            if (ii < 0)
                break;
            if (ii % 10000 == 0)
                cout << " --> Entry$= " << ii << endl;
            tree[s]->GetEntry(ii);
            if (s == 0)
                hid->Fill(dethit->event);
            //if (dethit->wave_max > 7 * dethit->ped_rms && dethit->q > 5)
            //if (dethit->wave_max > 3 * dethit->ped_rms && dethit->q > 5)
            //if (dethit->wave_max > 60 && dethit->q > 5)
            //if (dethit->wave_max > 150 && dethit->q > 6)
            //if (1)
            // x axis has a hit
            if (dethit->chn[0] != -999)
            {
                //if (dethit->wave_max > 10 * dethit->ped_rms && dethit->q > 10)
                if (dethit->wave_max > 7 * dethit->ped_rms && dethit->t < 300 && dethit->t > 270)
                {
                    hXHit[s]->Fill(dethit->chn[0]);
                    hXid[s]->Fill(dethit->event);

                    if (xid[s].empty())
                        xid[s].push_back(dethit->event);
                    else if (!IsMatched(dethit->event, xid[s]))
                        xid[s].push_back(dethit->event);
                }
            }

            // y axis has a hit
            if (dethit->chn[1] != -999)
            {
                //if (dethit->wave_max > 10 * dethit->ped_rms && dethit->q > 10)
                if (dethit->wave_max > 6 * dethit->ped_rms && dethit->t < 300 && dethit->t > 270)
                {
                    hYHit[s]->Fill(dethit->chn[1]);
                    hYid[s]->Fill(dethit->event);

                    if (yid[s].empty())
                        yid[s].push_back(dethit->event);
                    else if (!IsMatched(dethit->event, yid[s]))
                        yid[s].push_back(dethit->event);
                }
            }
        }
    }
    TH1I *hMMeff[3];
    TH1I *hlackMMeff[3];
    TH1I *hTrackereff[3];
    int matchedycounter;
    int matchedxcounter;
    int matchedandcounter;
    int matchedtotalcounter;
    vector<int> matchedyvec;
    vector<int> matchedxvec;
    vector<int> matchedandvec;
    vector<int> matchedorvec;
    for (int i = 0; i < 3; i++)
    {
        hMMeff[i] = new TH1I(Form("hMM%seff", XYname(i).Data()), "", 4, 0, 4);
        hlackMMeff[i] = new TH1I(Form("hlackMM%seff", XYname(i).Data()), "", 4, 0, 4);
        hTrackereff[i] = new TH1I(Form("hTracker%seff", XYname(i).Data()), "", 5, 0, 5);
        DrawMyHist(hMMeff[i], "chipID", "Counts", 1, 3);
        DrawMyHist(hlackMMeff[i], "chipID", "Counts", 1, 3);
        DrawMyHist(hTrackereff[i], "Num of firedchips", "Counts", 1, 3);
    }
    TH1I *hTrackerXYeff = new TH1I("hTrackerXYeff", "", 9, 0, 9);
    RecordStart = GetFirstTrig(hid);
    RecordEnd = GetLastTrig(hid);
    hid->GetXaxis()->SetRangeUser(RecordStart - 5, RecordEnd + 5);
    for (int s = RecordStart, jj = 0; s < RecordEnd + 1; s++)
    {
        DrawProcessbar(s - RecordStart, jj, RecordEnd + 1 - RecordStart);
        //if (s % 1000 == 0)
        //    cout << "The Entry No: " << s << endl;
        matchedxvec.clear();
        matchedyvec.clear();
        matchedandvec.clear();
        bool matchedx = 0;
        bool matchedy = 0;
        for (int i = 0; i < 4; i++)
        {
            matchedx = IsMatched(s, xid[i]);
            matchedy = IsMatched(s, yid[i]);
            if (matchedx)
            {

                matchedxvec.push_back(i);
                hMMeff[0]->Fill(i);
            }
            if (matchedy)
            {

                matchedyvec.push_back(i);
                hMMeff[1]->Fill(i);
            }
            if (matchedy && matchedx)
            {

                matchedandvec.push_back(i);
                hMMeff[2]->Fill(i);
            }
        }

        hTrackereff[0]->Fill(matchedxvec.size());
        if (matchedxvec.size())
        {
            matchedxcounter++;
            matchedtotalcounter++;
            if (matchedxvec.size() == 3)
                for (int j = 0; j < 4; j++)
                {

                    if (!IsMatched(j, matchedxvec))
                        hlackMMeff[0]->Fill(j);
                }
        }

        hTrackereff[1]->Fill(matchedyvec.size());
        if (matchedyvec.size())
        {
            matchedycounter++;
            matchedtotalcounter++;
            if (matchedyvec.size() == 3)
                for (int j = 0; j < 4; j++)
                {

                    if (!IsMatched(j, matchedyvec))
                        hlackMMeff[1]->Fill(j);
                }
        }

        hTrackereff[2]->Fill(matchedandvec.size());
        if (matchedandvec.size())
        {
            matchedandcounter++;
            if (matchedandvec.size() == 3)
                for (int j = 0; j < 4; j++)
                {

                    if (!IsMatched(j, matchedandvec))
                        hlackMMeff[2]->Fill(j);
                }
        }
        //if (matchedxvec.size() + matchedyvec.size() > 0)
        hTrackerXYeff->Fill(matchedxvec.size() + matchedyvec.size());
    }
    //NumofTriggers = RecordEnd - RecordStart + 1;

    NumofTriggers = RecordEnd - RecordStart + 1;
    int cnt = 0;
    for (int i = RecordStart; i < RecordEnd + 1; i++)
    {
        if (hid->GetBinContent(i) >= 1)
            cnt++;
    }
    RecordTriggers = cnt;
    int ValidTrackerXY = hTrackerXYeff->Integral(2, 100);
    int ValidTracker[3];
    ValidTracker[0] = hTrackereff[0]->Integral(2, 100);
    ValidTracker[1] = hTrackereff[1]->Integral(2, 100);
    ValidTracker[2] = hTrackereff[2]->Integral(2, 100);
    cout << "--> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl
         << "--> Trigger coincidence: " << NumofTriggers << endl
         << "--> Trigger Recorded by AGET: " << RecordTriggers << endl
         << " ** check ID lost please run function CheckTrackerID() **" << endl
         << "--> Valid tracker counter (or): " << ValidTrackerXY << ", GE: " << (double)ValidTrackerXY / RecordTriggers << endl
         << "--> Valid tracker x counter: " << ValidTracker[0] << ", GE: " << (double)ValidTracker[0] / RecordTriggers << endl
         << "--> Valid tracker y counter: " << ValidTracker[1] << ", GE: " << (double)ValidTracker[1] / RecordTriggers << endl
         << "--> Valid tracker (x&y) counter: " << ValidTracker[2] << ", GE: " << (double)ValidTracker[2] / RecordTriggers << endl;
    op << "--> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl
       << "--> Trigger coincidence: " << NumofTriggers << endl
       << "--> Trigger Recorded by AGET: " << RecordTriggers << endl
       << "--> Valid tracker counter (or): " << ValidTrackerXY << ", GE: " << (double)ValidTrackerXY / RecordTriggers << endl
       << "--> Valid tracker x counter: " << ValidTracker[0] << ", GE: " << (double)ValidTracker[0] / RecordTriggers << endl
       << "--> Valid tracker y counter: " << ValidTracker[1] << ", GE: " << (double)ValidTracker[1] / RecordTriggers << endl
       << "--> Valid tracker (x&y) counter: " << ValidTracker[2] << ", GE: " << (double)ValidTracker[2] / RecordTriggers << endl;
    TLatex *la;
    TCanvas *c;
    int ccounter = 0;

    TString preName = fAnaName.Copy();
    preName = preName.Remove(preName.Length() - 5, 5);

    TString PngName;

    //** draw ID distribution
    //*
    c = new TCanvas(Form("c%d", ccounter), Form("c%d", ccounter), 1500, 1200);
    c->Divide(4, 4);
    for (int i = 0; i < 4; i++)
    {
        c->cd(i + 1);
        hXHit[i]->Draw();
        c->cd(i + 5);
        hYHit[i]->Draw();
        c->cd(i + 9);
        hXid[i]->Draw();
        c->cd(i + 13);
        hYid[i]->Draw();
    }
    PngName = preName + "TrackerIDdist.png";
    c->SaveAs(PngName);
    ccounter++;

    //** draw tracker XYeff
    //*
    c = cdC(ccounter++);
    int noevent = RecordTriggers - ValidTrackerXY;
    hTrackerXYeff->SetBinContent(hTrackerXYeff->FindBin(0), noevent);
    DrawMyHist(hTrackerXYeff, "Num of firedchips", "Counts", 1, 3);
    hTrackerXYeff->Draw();
    c->Update();
    SetEffstats(c, hTrackerXYeff, NumofTriggers, RecordTriggers, ValidTrackerXY);

    double Trackereff[8 + 1] = {0};

    Trackereff[0] = hTrackerXYeff->GetBinContent(1) / RecordTriggers;
    sprintf(buff, "No chip response,Ratio=%.1f%%", Trackereff[0] * 100);
    la = DrawMyLatex(buff, 0.25, 0.7, 42, 0.05, 1);
    la->Draw("same");
    for (int i = 1; i <= 8; i++)
    {

        Trackereff[i] = hTrackerXYeff->GetBinContent(i + 1) / RecordTriggers;
        sprintf(buff, "Nchips=%d,Ratio=%.1f%%,#color[2]{(valid)%.1f%%}", i, Trackereff[i] * 100, Trackereff[i] * 100 * RecordTriggers / ValidTrackerXY);
        la = DrawMyLatex(buff, 0.25, 0.7 - 0.05 * i, 42, 0.04);
        la->Draw("same");
    }
    PngName = preName + "TrackerXYeff.png";
    c->SaveAs(PngName);

    for (int s = 0; s < 3; s++)
    {
        //** draw tracker eff
        //*
        c = cdC(ccounter++);
        c->Clear();
        int noevent = RecordTriggers - ValidTracker[s];
        hTrackereff[s]->SetBinContent(hTrackereff[s]->FindBin(0), noevent);

        DrawMyHist(hTrackereff[s], "Num of MM", "Counts", 1, 3);
        hTrackereff[s]->Draw();
        c->Update();
        SetEffstats(c, hTrackereff[s], NumofTriggers, RecordTriggers, ValidTracker[s]);
        double Trackereff[4 + 1] = {0};

        Trackereff[0] = hTrackereff[s]->GetBinContent(1) / RecordTriggers;
        sprintf(buff, "No MM response,Ratio=%.1f%%", Trackereff[0] * 100);
        la = DrawMyLatex(buff, 0.25, 0.7, 42, 0.05, 1);
        la->Draw("same");
        for (int i = 1; i <= 4; i++)
        {

            Trackereff[i] = hTrackereff[s]->GetBinContent(i + 1) / RecordTriggers;
            sprintf(buff, "NMM=%d,Ratio=%.1f%%,#color[2]{(valid)%.1f%%}", i, Trackereff[i] * 100, Trackereff[i] * 100 * RecordTriggers / ValidTracker[s]);
            //sprintf(buff, "Num of MM=%d,Ratio=%.1f%%", i + 1, Trackereff[i] * 100);
            la = DrawMyLatex(buff, 0.25, 0.7 - 0.07 * i, 42, 0.05);
            //la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
            la->Draw("same");
        }
        PngName = preName + Form("Tracker%seff.png", XYname(s).Data());
        c->SaveAs(PngName);

        //** draw MM eff
        //*
        c = cdC(ccounter++);
        DrawMyHist(hMMeff[s], "ID of MM", "Counts", 1, 3);
        hMMeff[s]->Draw();
        c->Update();
        SetEffstats(c, hMMeff[s], NumofTriggers, RecordTriggers, ValidTracker[s]);

        double MMeff[4];
        for (int i = 0; i < 4; i++)
        {

            MMeff[i] = hMMeff[s]->GetBinContent(i + 1) / ValidTracker[s];
            sprintf(buff, "MMID=%d,Ratio=%.1f%%", i, MMeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.7 - 0.07 * i, 42, 0.05, 2);
            la->Draw("same");
        }
        PngName = preName + Form("MM%seff.png", XYname(s).Data());
        c->SaveAs(PngName);

        //** draw MM true eff
        //*
        int allfires = allfires = hTrackereff[s]->GetBinContent(4 + 1);
        c = cdC(ccounter++);
        DrawMyHist(hlackMMeff[s], "ID of MM", "Counts", 1, 3);
        hlackMMeff[s]->Draw();
        c->Update();
        SetEffstats(c, hlackMMeff[s], NumofTriggers, RecordTriggers, ValidTracker[s]);
        double MMtrueeff[4];
        for (int i = 0; i < 4; i++)
        {

            MMtrueeff[i] = allfires / (allfires + hlackMMeff[s]->GetBinContent(i + 1));
            sprintf(buff, "MMID=%d, Eff=%.1f%%", i, MMtrueeff[i] * 100);
            la = DrawMyLatex(buff, 0.3, 0.7 - 0.07 * i, 42, 0.05, 2);
            la->Draw("same");
        }
        PngName = preName + Form("MM%strueeff.png", XYname(s).Data());
        c->SaveAs(PngName);
    }
}

void DefineRBCRhist()
{
    htheta = new TH1F("theta", "theta", 100, 0, 0.5);
    hphi = new TH1F("phi", "phi", 200, -4, 4);
    hT0map = new TH2F("T0map", "", 100, -200, 200, 100, -200, 200);
    hFTOFmap = new TH2F("FTOFmap", "", 100, -200, 200, 100, -200, 200);
}

void DrawRBCRresults(TString cname)
{

    TCanvas *cc = cdC(1, 1800, 600);
    cc->Clear();
    cc->Divide(4, 1);
    cc->cd(1);
    htheta->Draw();
    cc->cd(2);
    hphi->Draw();
    cc->cd(3);
    hT0map->Draw("colz");
    cc->cd(4);
    hFTOFmap->Draw("colz");
    cc->SaveAs(cname.Data());
}
void Definedecodehist()
{
    //--------------------------
    // 定义histogram & graph
    int nbin;

    //事例填图
    for (int i = 0; i < NumMM; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            /*
            */
            nbin = 64;
            fAGETch[i][j] = new TH1F(Form("fAGETch%d_%d", i, j), Form("ID of fired AGET channel of %s for %s", detName[i], xyName[j]), 100, 0, nbin);
            fCharge[i][j] = new TH1F(Form("fQ%d_%d", i, j), Form("AGET channel charge of %s for %s", detName[i], xyName[j]), 500, 0, 500);
            fAGETchNum[i][j] = new TH1F(Form("fAGETchNum%d_%d", i, j), Form("number of fired AGET channel of %s for %s", detName[i], xyName[j]), 30, 0, 30);
            for (int k = 0; k < 2; k++)
            {

                fChargesum[k][i][j] = new TH1F(Form("fQsum%d_%d", i, j), Form("cluster charge of %s for %s cut%d", detName[i], xyName[j], k), 1e3, 0, 1000);

                fClustersize[k][i][j] = new TH1F(Form("fclustersize%d_%d", i, j), Form("Cluster size of %s for %s cut%d", detName[i], xyName[j], k), 30, 0, 30);
                nbin = 384;
                fClustermean[k][i][j] = new TH1F(Form("fclustermean%d_%d", i, j), Form("Cluster mean of %s for %s cut%d", detName[i], xyName[j], k), nbin, 0, nbin);
                fClusterNum[k][i][j] = new TH1F(Form("fclusternum%d_%d", i, j), Form("Number of likely cases of cluster of %s for %s cut", detName[i], xyName[j], k), 100, 0, 100);
            }

            fPeakingTime[i][j] = new TH1F(Form("fT%d", i), Form("Peaking time of %s", detName[i]), 512, 0, 512);
            fQT[i][j] = new TH2F(Form("fQT%d", i), Form("Q vs. T of %s", detName[i]), 100, 0, 500, 100, 0, 512);
            fQT[i][j]->GetXaxis()->SetTitle("Q");
            fQT[i][j]->GetYaxis()->SetTitle("T");
        }
    }
}
void drawdecoderesults(TString cname, int id)
{

    TCanvas *cc = new TCanvas(cname);
    cc->SetTitle(detName[id]);
    cc->Clear();
    cc->Divide(4, 4);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            cc->cd(1 + 4 * j + 8 * i);
            fChargesum[i][id][j]->Draw("");
            cc->cd(2 + 4 * j + 8 * i);
            fClustersize[i][id][j]->Draw("");
            cc->cd(3 + 4 * j + 8 * i);
            fClustermean[i][id][j]->Draw();
            cc->cd(4 + 4 * j + 8 * i);
            fClusterNum[i][id][j]->Draw();
        }
    }
    /*
    cc->cd(1);
    fAGETch[id][0]->Draw("colz");
    cc->cd(2);
    fCharge[id][0]->Draw("");
    cc->cd(3);
    fChargesum[id][0]->Draw("");
    cc->cd(4);
    fAGETchNum[id][0]->Draw("");
    cc->cd(5);
    fAGETch[id][1]->Draw("colz");
    cc->cd(6);
    fCharge[id][1]->Draw("");
    cc->cd(7);
    fChargesum[id][1]->Draw("");
    cc->cd(8);
    fAGETchNum[id][1]->Draw("");
    cc->cd(9);
    fClustersize[id][0]->Draw("");
    cc->cd(10);
    fClustermean[id][0]->Draw();
    //fPeakingTime[id]->Draw("");
    cc->cd(11);
    fClusterNum[id][0]->Draw();
    //fCharge[id]->Draw("");
    cc->cd(12);
    //fQT[id]->Draw("colz");
    cc->cd(13);
    fClustersize[id][1]->Draw("");
    cc->cd(14);
    fClustermean[id][1]->Draw();
    cc->cd(15);
    fClusterNum[id][1]->Draw();
    */
}
#if 0
void DecodeAGET(TString fAnaName, TString fDecName, int force)
{
    TString fpath = GetPath(fAnaName);
    ofstream op[2];
    op[0].open(Form("%s/Tracker_xdecode.dat", fpath.Data()));
    op[1].open(Form("%s/Tracker_ydecode.dat", fpath.Data()));

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fDecName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fDecName << " is exist!" << endl;
            return;
        }
    }

    TFile *decFile = new TFile(fDecName, "recreate");
    if (!decFile->IsOpen())
    {
        cout << "ERROR! " << fDecName << " cant open!" << endl;
        return;
    }
    Definedecodehist();
    TTree *dtree[NumMM];
    MMHitData hittemp;
    MMHit MMhit;
    for (int s = 0; s < NumMM; s++)
    {
        dtree[s] = new TTree(Form("tree%d", s), Form("tree%d", s));
        dtree[s]->Branch("MMhit", &MMhit);
    }
    hit_data_dec_total *trackerhit = new hit_data_dec_total("asic2det_64to384.csv");
    //return;
    AGETHit *dethit;
    dethit = new AGETHit;

    vector<int> chn[2];
    vector<double> q[2];

    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[NumMM];
    int allmin = 9999999;
    int allmax = 0;
    int min[NumMM];
    int max[NumMM];
    for (int s = 0; s < NumMM; s++)
    {
        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        tree[s]->Print();
        tree[s]->SetBranchAddress("dethit", &dethit);
        int nentries = tree[s]->GetEntriesFast();
        cout << " --> nentries= " << nentries << endl;
        tree[s]->GetEntry(0);
        min[s] = dethit->event;
        tree[s]->GetEntry(nentries - 1);
        max[s] = dethit->event;
        cout << "--> "
             << "Tracker" << s << ", Events: " << nentries << ", tigger range: " << min[s] << "to " << max[s] << endl;
    }
    allmin = TMath::MinElement(NumMM, min);
    allmax = TMath::MaxElement(NumMM, max);
    cout << "--> Searching trigger from " << allmin << " to " << allmax << endl;
    int iEntry[NumMM] = {0};
    vector<AGETHit> AGET_chn[2];
    vector<double> validclusterm;
    int xymatched[NumMM][2] = {0};
    for (int trigid = allmin; trigid < allmax && trigid >= 0 && trigid < 100e4; trigid++)
    {
        if (trigid % 100 == 0)
            cout << "--> Searching trigger id = " << trigid << endl;

        for (int s = 0; s < NumMM; s++)
        {
            MMhit.initial();
            int nrange = 200;
            int nentries = tree[s]->GetEntriesFast();
            int start = (iEntry[s] - nrange < 0) ? 0 : iEntry[s] - nrange;
            int end = (iEntry[s] + nrange > nentries) ? nentries : iEntry[s] + nrange;

            vector<int>().swap(chn[0]);
            vector<int>().swap(chn[1]);
            vector<double>().swap(q[0]);
            vector<double>().swap(q[1]);
            vector<AGETHit>().swap(AGET_chn[1]);
            vector<AGETHit>().swap(AGET_chn[1]);
            //cout<<"MM"<<detName[s]<<", iEntry: "<<iEntry[s]<<", start: "<<start<<", end: "<<end<<endl;
            for (int i = start; i < end; i++)
            {
                tree[s]->GetEntry(i);
                if (dethit->event == trigid)
                {
                    iEntry[s] = i;
                    if (dethit->wave_max > 150 && dethit->q > 6)
                    //if (dethit->wave_max > 7 * dethit->ped_rms && dethit->q > 5)
                    {
                        for (int j = 0; j < 2; j++)
                        {
                            if (dethit->chn[j] >= 0 && dethit->chn[j] < 64)
                            {
                                chn[j].push_back(dethit->chn[j]);
                                q[j].push_back(dethit->q);
                                AGET_chn[j].push_back(*dethit);
                                fAGETch[s][j]->Fill(chn[j].back());
                                fCharge[s][j]->Fill(q[j].back());
                            }
                        }
                    }
                }
            }
            int xyvalid = 0;
            for (int j = 0; j < 2; j++)
            {
                if (!chn[j].size() || !q[j].size())
                    continue;
                op[j] << xyName[j] << ": " << endl;
                op[j] << "AGET CHn: ";
                for (int k = 0; k < chn[j].size(); k++)
                {
                    op[j] << " " << chn[j][k];
                }
                op[j] << endl;

                op[j] << "AGET Charge: ";
                for (int k = 0; k < q[j].size(); k++)
                {
                    op[j] << " " << q[j][k];
                }
                op[j] << endl;
                xymatched[s][j]++;
                fAGETchNum[s][j]->Fill(chn[j].size());
                hittemp.initial();
                trackerhit->clear();
                if (trackerhit->hit_dec(chn[j], q[j]))
                {

                    trackerhit->result_sort();
                    validclusterm.clear();
                    //cout << xyName[j] << " hit pos rebuild success!" << endl;

                    op[j] << "Decode Charge: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_amp.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode clustermean: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_strip.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode AGET_chn_num_use: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->aget_chn_num_use.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode cluster size: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_strip_num.at(k);
                    }
                    op[j] << endl;

                    //for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    for (int k = 0; k < 1; k++)
                    {
                        //hittemp.chnum = xchn;
                        //hittemp.chq = xq;
                        hittemp.AGET_chn = AGET_chn[0];
                        hittemp.q = trackerhit->hit_amp.at(k);
                        hittemp.stripnum = trackerhit->hit_strip_num.at(k);
                        hittemp.AGET_chn_num_use = trackerhit->aget_chn_num_use.at(k);
                        hittemp.clustermean[j] = trackerhit->hit_strip.at(k);
                        MMhit.MMchn[j].push_back(hittemp);
                        //cout << "hit_strip=" << trackerhit->hit_strip.at(k) << endl
                        //     << "hit_charge=" << trackerhit->hit_amp.at(k) << endl
                        //     << "hit_strip_num=" << trackerhit->hit_strip_num.at(k) << endl
                        //     << "AGET_chn_num_use=" << trackerhit->aget_chn_num_use.at(k) << endl
                        //     << endl;
                        //if(IsTrackClusterEffective(hittemp))
                        fChargesum[0][s][j]->Fill(hittemp.q);
                        fClustersize[0][s][j]->Fill(hittemp.stripnum);
                        fClustermean[0][s][j]->Fill(hittemp.clustermean[j]);
                        if (IsTrackClusterEffective(hittemp))
                        {

                            validclusterm.push_back(hittemp.clustermean[j]);
                            fChargesum[1][s][j]->Fill(hittemp.q);
                            fClustersize[1][s][j]->Fill(hittemp.stripnum);
                            fClustermean[1][s][j]->Fill(hittemp.clustermean[j]);
                        }
                    }

                    fClusterNum[0][s][j]->Fill(trackerhit->hit_strip.size());
                    if (validclusterm.size())
                    {

                        fClusterNum[1][s][j]->Fill(validclusterm.size());
                        //fClusterNum[s][j]->Fill(trackerhit->hit_strip.size());
                        MMhit.xyflag = j;
                        xyvalid++;
                    }
                }
                else
                    op[j] << "Decode failed" << endl;
            }
            if (xyvalid > 0)
            {

                MMhit.MMID = s;
                MMhit.event = trigid;
                if (xyvalid >= 2)
                    MMhit.xyflag = xyvalid;
                dtree[s]->Fill();
            }
        }
    }
    /*
    for (int s = 0; s < NumMM; s++)
    {
        cout << "--> "
             << "Tracker" << s << ", tigger range: " << min[s] << "to " << max[s] << endl
             << " matched x:" << xymatched[s][0] << ", matched y:" << xymatched[s][1] << endl;
        drawdecoderesults(Form("c%d", s), s);
    }
*/
    for (int i = 0; i < NumMM; i++)
    {
        decFile->WriteTObject(dtree[i]);
    }
    decFile->Close();
    anaFile->Close();
    cout << "--> Decode results has been stored to " << fDecName << endl
         << endl;
}
#endif
void DecodeAGET4Align(TString fAnaName, TString fDecName, int force)
{
    TString fpath = GetPath(fAnaName);
    ofstream op[2];
    ofstream oplog;

    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fDecName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fDecName << " is exist!" << endl;
            return;
        }
    }

    TFile *decFile = new TFile(fDecName, "recreate");
    if (!decFile->IsOpen())
    {
        cout << "ERROR! " << fDecName << " cant open!" << endl;
        return;
    }
    int event;
    PositionData detector_data[NumMM];
    TTree *dtree = new TTree("tree", "tree");
    dtree->Branch("event", &event, "event/I");
    for (int s = 0; s < NumMM; s++)
    {
        detector_data[s] = PositionData();
        string dec_name = "dector" + to_string(s);
        dtree->Branch(dec_name.c_str(), &detector_data[s].sig, "sig/L:sig_x/L:sig_y/L:x/D:y/D:z/D:hit_strip_x/L:hit_strip_y/L:hit_amp_x/D:hit_amp_y/D:hit_chn_x/L:hit_chn_y/L:x_nhits/L:x_other[5]/D:y_nhits/L:y_other[5]/D");
    }
    //dtree->Print();
    hit_data_dec_total *trackerhit = new hit_data_dec_total("asic2det_64to384.csv");
    //return;
    AGETHit *dethit;
    dethit = new AGETHit;
    vector<AGETHit> AGETData[NumMM];

    vector<int> chn[2];
    vector<double> q[2];

    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[NumMM];
    int allmin = 9999999;
    int allmax = 0;
    int min[NumMM];
    int max[NumMM];
    for (int s = 0; s < NumMM; s++)
    {
        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        //tree[s]->Print();
        tree[s]->SetBranchAddress("dethit", &dethit);
        int nentries = tree[s]->GetEntriesFast();
        cout << " --> nentries= " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {

            tree[s]->GetEntry(i);
            AGETData[s].push_back(*dethit);
            //cout<<"AGETData[s].back().event: "<<AGETData[s].back().event<<endl;
        }
        tree[s]->GetEntry(0);
        min[s] = dethit->event;
        tree[s]->GetEntry(nentries - 1);
        max[s] = dethit->event;
        cout << "--> "
             << "Tracker" << s << ", Events: " << nentries << ", tigger range: " << min[s] << "to " << max[s] << endl;
    }
    allmin = TMath::MinElement(NumMM, min);
    allmax = TMath::MaxElement(NumMM, max);
    cout << "--> Searching trigger from " << allmin << " to " << allmax << endl;
    int iEntry[NumMM] = {0};
    op[0].open(Form("%s/Tracker_xdecode.dat", fpath.Data()));
    op[1].open(Form("%s/Tracker_ydecode.dat", fpath.Data()));
    oplog.open(Form("%s/Trackerdecode_log.dat", fpath.Data()));
    cout << "Create logfile: " << Form("%s/Trackerdecode_log.dat", fpath.Data()) << endl;
    //for (int trigid = 0; trigid < allmax && trigid >= 0 && trigid < 100e4; trigid++)
    for (int trigid = allmin, pb = 0; trigid < allmax && trigid >= 0 && trigid < 100e4; trigid++)
    {
        DrawProcessbar(trigid - allmin, pb, allmax - allmin);
        //if (trigid % 100 == 0)
        //    cout << "--> Searching trigger id = " << trigid << endl;
        int validMM = 0;
        event = -999;
        for (int s = 0; s < NumMM; s++)
        //for (int s = 2; s < 3; s++)
        {
            //oplog<<"progress check?" <<endl;
            detector_data[s].Initial();
            int nrange = 200;
            int nentries = tree[s]->GetEntriesFast();
            int start = (iEntry[s] - nrange < 0) ? 0 : iEntry[s] - nrange;
            int end = (iEntry[s] + nrange > nentries) ? nentries : iEntry[s] + nrange;
            //int start = 0;
            //int end = nentries;

            vector<int>().swap(chn[0]);
            vector<int>().swap(chn[1]);
            vector<double>().swap(q[0]);
            vector<double>().swap(q[1]);
            //cout<<"MM"<<detName[s]<<", iEntry: "<<iEntry[s]<<", start: "<<start<<", end: "<<end<<endl;
            if (s == 0 && trigid > 0 && trigid < 1e3)
                oplog << "MM" << detName[s] << ", iEntry: " << iEntry[s] << ", start: " << start << ", end: " << end << endl;
            for (int i = start; i < end; i++)
            {
                //tree[s]->GetEntry(i);
                //cout<<"AGETData[s][i].event, trigid"<<AGETData[s][i].event<<","<<trigid<<endl;
                //if (s == 0 && trigid > 0 && trigid < 1e3)
                //if (s == 0 && trigid > 79e3 && trigid < 81e3)
                if (AGETData[s][i].event == trigid)
                {
                    iEntry[s] = i;

                    //if (AGETData[s][i].wave_max > 7 * AGETData[s][i].ped_rms && AGETData[s][i].q > 5)
                    //{
                    for (int j = 0; j < 2; j++)
                    {
                        switch (j + 1)
                        {
                        case 1:
                            if (
                                AGETData[s][i].wave_max > 7 * AGETData[s][i].ped_rms
                                //AGETData[s][i].wave_max > 10 * AGETData[s][i].ped_rms && AGETData[s][i].q > 10
                                && AGETData[s][i].t > 270 && AGETData[s][i].t < 290 && AGETData[s][i].chn[j] >= 0 && AGETData[s][i].chn[j] < 64)
                            {
                                chn[j].push_back(AGETData[s][i].chn[j]);
                                q[j].push_back(AGETData[s][i].wave_max);
                            }
                            break;
                        case 2:
                            if (
                                AGETData[s][i].wave_max > 6 * AGETData[s][i].ped_rms
                                //AGETData[s][i].wave_max > 10 * AGETData[s][i].ped_rms && AGETData[s][i].q > 10
                                && AGETData[s][i].t > 270 && AGETData[s][i].t < 290 && AGETData[s][i].chn[j] >= 0 && AGETData[s][i].chn[j] < 64)
                            {
                                chn[j].push_back(AGETData[s][i].chn[j]);
                                q[j].push_back(AGETData[s][i].wave_max);
                                //cout<<"chn, q: " <<chn[j].back()<<","<<q[j].back()<<endl;
                            }
                            break;

                        default:
                            break;
                        }
                    }
                    //end = i+200;
                    //}
                }
            }
            for (int j = 0; j < 2; j++)
            {
                if (!chn[j].size() || !q[j].size())
                    continue;
                op[j] << xyName[j] << ": " << endl;
                op[j] << "AGET CHn: ";
                for (int k = 0; k < chn[j].size(); k++)
                {
                    op[j] << " " << chn[j][k];
                }
                op[j] << endl;

                op[j] << "AGET Charge: ";
                for (int k = 0; k < q[j].size(); k++)
                {
                    op[j] << " " << q[j][k];
                }
                op[j] << endl;

                trackerhit->clear();
                if (trackerhit->hit_dec(chn[j], q[j]))
                {

                    validMM++;
                    trackerhit->result_sort();
                    //cout << xyName[j] << " hit pos rebuild success!" << endl;

                    op[j] << "Decode Charge: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_amp.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode clustermean: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_strip.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode AGET_chn_num_use: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->aget_chn_num_use.at(k);
                    }
                    op[j] << endl;

                    op[j] << "Decode cluster size: ";
                    for (int k = 0; k < trackerhit->hit_strip.size(); k++)
                    {
                        op[j] << " " << trackerhit->hit_strip_num.at(k);
                    }
                    op[j] << endl;
                    if (j == 0)
                    {
                        detector_data[s].x = trackerhit->hit_strip.at(0);
                        detector_data[s].hit_chn_num_x = trackerhit->aget_chn_num_use.at(0);
                        detector_data[s].hit_strip_num_x = trackerhit->hit_strip_num.at(0);
                        detector_data[s].hit_amp_x = trackerhit->hit_amp.at(0);
                        detector_data[s].sig_x = 1;
                        int cnt = 0;
                        for (int k = 1; k < trackerhit->hit_strip.size() && cnt < 5; k++)
                        {
                            detector_data[s].x_other[cnt] = trackerhit->hit_strip.at(k);
                            cnt++;
                        }
                        //TODO
                        detector_data[s].x_nhits = cnt;
                    }
                    if (j == 1)
                    {
                        detector_data[s].y = trackerhit->hit_strip.at(0);
                        detector_data[s].hit_chn_num_y = trackerhit->aget_chn_num_use.at(0);
                        detector_data[s].hit_strip_num_y = trackerhit->hit_strip_num.at(0);
                        detector_data[s].hit_amp_y = trackerhit->hit_amp.at(0);
                        detector_data[s].sig_y = 1;
                        int cnt = 0;
                        for (int k = 1; k < trackerhit->hit_strip.size() && cnt < 5; k++)
                        {
                            detector_data[s].y_other[cnt] = trackerhit->hit_strip.at(k);
                            cnt++;
                        }
                        //TODO
                        detector_data[s].y_nhits = cnt;
                    }
                    detector_data[s].z = Zpos[s];
                }
                else
                    op[j] << "Decode failed" << endl;
            }

            if (detector_data[s].sig_y && detector_data[s].sig_x)
                detector_data[s].sig = 1;
        }
        if (validMM > 0)
        {

            event = trigid;
            dtree->Fill();
        }
    }
    /*
    for (int s = 0; s < NumMM; s++)
    {
        cout << "--> "
             << "Tracker" << s << ", tigger range: " << min[s] << "to " << max[s] << endl
             << " matched x:" << xymatched[s][0] << ", matched y:" << xymatched[s][1] << endl;
        drawdecoderesults(Form("c%d", s), s);
    }
*/
    //oplog.close();
    //op[0].close();
    //op[1].close();
    decFile->WriteTObject(dtree);
    decFile->Close();
    anaFile->Close();
    cout << "--> Decode results has been stored to " << fDecName << endl
         << endl;
}
bool RebuildCRAngle(vector<double> x, vector<double> y, vector<double> z, RBCRdata &RBCRpos, double fit_R2_cut = 0.9)
{
    int N = z.size();
    if (N < 3)
        return false;
    //* 0-x, 1-y, z=0;
    double expT0[2] = {0};
    double expFTOF[2] = {0};
    double delta[3] = {0};
    vector<double> xy[2];
    xy[0] = x;
    xy[1] = y;
    double par[2][2];
    TVector3 v1;
    for (int i = 0; i < 2; i++)
    {
        double fit_R2 = linearfit(z, xy[i], par[i]);
        cout << " correlation index: " << fit_R2 << endl;
        if (fit_R2 > fit_R2_cut)
        {

            expT0[i] = par[i][0] * Zpos[T0id] + par[i][1];
            expFTOF[i] = par[i][0] * Zpos[FTOFid] + par[i][1];
            delta[i] = par[i][0] * Zpos[MM3id] + par[i][1] - (par[i][0] * Zpos[MM0id] + par[i][1]);
        }
    }
    delta[2] = Zpos[MM3id] - Zpos[MM0id];
    v1.SetX(-1 * delta[0]);
    v1.SetY(-1 * delta[1]);
    v1.SetZ(-1 * delta[2]);
    RBCRpos.initial();
    if (expT0[0] && expT0[1])
    {
        RBCRpos.T0detX = expT0[0];
        RBCRpos.T0detY = expT0[1];
        RBCRpos.FTOFdetX = expFTOF[0];
        RBCRpos.FTOFdetY = expFTOF[1];
        RBCRpos.Mupx = -1 * delta[0] / v1.Mag();
        RBCRpos.Mupy = -1 * delta[1] / v1.Mag();
        RBCRpos.Mupz = -1 * delta[2] / v1.Mag();
        RBCRpos.Mutheta = v1.Theta();
        RBCRpos.Muphi = v1.Phi();
        cout << "theta: phi" << v1.Theta() << "\t" << v1.Phi() << endl;
        return true;
    }
    return false;
}
double RebuildT0Track(double np, TVector3 Inpos, TVector3 Indir, int ID = 1)
{
    double Inx = Inpos.X();
    double Iny = Inpos.Y();
    double Inz = Inpos.Z();
    double InPx = Indir.X();
    double InPy = Indir.Y();
    double InPz = Indir.Z();
    swap(Inx, Inz);
    swap(InPx, InPz);
    //swap(Iny,Inz);
    //swap(InPy,InPz);

    /*
           double Ox, Oy, Oz;
           Oz = A1z;
           Ox = (A1z - Inz) / pz * px + Inx;
        //cout<<Ox<<endl;
        Oy = (A1z - Inz) / pz * py + Iny;
        return sqrt(Ox * Ox + Oy * Oy + Oz * Oz);
        //cout << "A1z: " << A1z<<"\t"<<Inz<<"\t"<<pz<<"\t"<<px<<"\t"<<Inx << endl;
        */
    double detpos = Zpos[T0];
    double Ax0, Ay0, Az0;
    double px, py, pz;
    //double swap;
    if (ID == 1 || ID == 3)
    {
        Az0 = -1 * (ID - 2) * 90;
        Ax0 = detpos;
        Ay0 = 0;
        px = InPx;
        py = InPy;
        pz = -1 * (ID - 2) * InPz;
    }
    if (ID == 2 || ID == 4)
    {
        Az0 = -1 * (ID - 3) * 90;
        Ax0 = detpos;
        Ay0 = 0;
        px = InPx;
        py = InPz;
        pz = -1 * (ID - 3) * InPy;
        swap(Iny, Inz);
        //swap = Iny;
        //Iny = Inz;
        //Inz = swap;
    }
    double thetaC = TMath::ACos(1 / np);
    double theta = TMath::ACos(pz / TMath::Sqrt(px * px + py * py + pz * pz)); //angle between momentum direction of mu and z axis;

    double Ax, Ay, Az;
    Az = Az0;
    double dx, dy, dz;
    dz = Az0 - Inz;
    //return TMath::Abs(dz/TMath::Cos(theta - thetaC));

    double Rdxdy;
    Rdxdy = dz * TMath::Tan(theta - thetaC);
    dx = Rdxdy / sqrt(px * px + py * py) * px;
    dy = Rdxdy / sqrt(px * px + py * py) * py;
    //Ax = Inx + dx - detpos;
    Ax = dx;
    Ay = Iny + dy;
    if (ID == 2 || ID == 4)
        swap(Ay, Az);
    //cout << "Input: " <<ID<<"\t"<< Inx << "\t" << Iny << "\t" << Inz << "\t" << InPx << "\t" << InPy << "\t" << InPz << endl;
    //cout<<"Hit pos="<<Ax<<", "<<Ay<<", "<<Az<<endl;
    //cout<<"Hit pos="<<Ax<<", "<<Ay<<", "<<Az<<endl;
    return sqrt(dx * dx + dy * dy + dz * dz);
    //return sqrt(Ax * Ax + Ay * Ay + Az * Az);

    //return -1;
}
vector<double> TAcorrection(TString name, int iter, vector<double> TT, vector<double> AA, int rbt = 10, int rbA = 12, double initialtL = -2e3, double initialtR = 3e3, double initialAL = -300, double initialAR = 300)
{
    bint = (initialtR - initialtL) / 0.5;
    cout << "bint: " << bint << endl;
    cout << "tL, tR: " << initialtL << "\t" << initialtR << endl;
    ht = new TH1D("ht", ";Time (ns);Counts", bint, initialtL, initialtR);
    htcor = new TH1D("htcor", ";Time (ns);Counts", bint, initialtL, initialtR);
    hL = new TH1D("hL", ";reTrack (mm);Counts", 3e3, initialAL, initialAR);
    hLT = new TH2D("hLT", ";reTrack;TR (ns)", 3e3, initialAL, initialAR, bint, initialtL, initialtR);

    TCanvas *cc;
    int CNum = 0;
    TH1F *htfit;
    TH1F *hAfit;
    TF1 *fitT;
    TF1 *fitA;
    TF1 *fitAT = new TF1("fitAT", "0", initialAL, initialAR);

    double Treserve[200000];
    double TAcor = 0;
    double Amean = 0;
    double Asigma = 0;
    double Tmean = 0;
    double Tsigma = 0;
    double fTL = 0;
    double fTR = 0;
    vector<double> Tcor;
    cc = cdC(CNum++);
    for (int s = 0; s < iter; s++)
    {
        if (s != 0)
        {
            //if (Amean - 3 * Asigma <= 0)
            //    fitAT = profilefit(hLT, rbA, rbt * 4, Tmean - 6 * Tsigma, Tmean + 6 * Tsigma, 0.1, Amean + 3 * Asigma, Form("%s/%d", path.Data(), s-1));
            //else
            fitAT = profilefit(hLT, rbA, rbt * 4, Tmean - 6 * Tsigma, Tmean + 6 * Tsigma, Amean - 3 * Asigma, Amean + 3 * Asigma, Form("%s_%d", name.Data(), s - 1));
            if (!fitAT)
            {
                cout << " the profilefit is failed! " << endl;
                return Tcor;
            }
            ht->Reset();
            hLT->Reset();
        }
        for (int i = 0; i < TT.size(); i++)
        {
            if (s == 0)
            {
                //cout << "AA & TT: " << AA[i] << "\t" << TT[i] << endl;
                hL->Fill(AA[i]);
                Treserve[i] = TT[i];
                fTL = initialtL;
                fTR = initialtR;
            }
            else
            {

                TAcor = fitAT->Eval(AA[i]);
                //cout << "TAcor: " << TAcor << endl;
                Treserve[i] = Treserve[i] - TAcor;
                fTL = -1e3;
                fTR = 1e3;
            }
            ht->Fill(Treserve[i]);
            hLT->Fill(AA[i], Treserve[i]);
        }

        // * draw Tracklength
        cc = cdC(100);
        if (s == 0)
        {
            //DrawMyHist(hL, "", "", 1, 3);
            hL->Draw();
            fitA = gausfit(hL, 100, 3, 3, rbA * 4, initialAL, initialAR);
            //return;
            if (!fitA)
            {
                cout << "the hAfit hist is NULL " << endl;
                return Tcor;
            }
            //fitA = (TF1 *)hAfit->GetFunction("fitU");
            Amean = fitA->GetParameter(1);
            Asigma = fitA->GetParameter(2);
            cout << "Amean=" << Amean << ",\tAsigma=" << Asigma << endl;
            cc->SaveAs(Form("%s.png", name.Data()));
        }

        // * draw T
        cc = cdC(101);
        cc->Clear();
        //DrawMyHist(ht, "", "", 1, 3);
        ht->Draw();
        fitT = gausfit(ht, 100, 3, 3, rbt * 4, fTL, fTR);
        //fitT = gausfit(ht, 0.1, 3, 3, rbt * 20, fTL, fTR);
        if (!fitT)
        {
            cout << "the htfit hist is NULL " << endl;
            return Tcor;
        }
        //fitT = (TF1 *)htfit->GetFunction("fitU");
        Tmean = fitT->GetParameter(1);
        Tsigma = fitT->GetParameter(2);

        ht->SetNdivisions(505);
        sprintf(buff, "#sigma=%.0fps", Tsigma);
        TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
        l->Draw();
        cout << "Tmean=" << Tmean << ",\tTsigma=" << Tsigma << endl;
        cc->SaveAs(Form("%s_TR_cor%d.png", name.Data(), s));

        //* draw T vs track
        cc = cdC(102);
        cc->Clear();
        //DrawMy2dHist(hLT, "", "", 1, 2);
        //ht->Rebin(2);
        hLT->Draw("colz");
        cc->Modified();
        cc->Update();

        //fit = gausfit(ht, 20e-3, 3, 3, 1, fTL, fTR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        cc->SaveAs(Form("%s_hAT%d.png", name.Data(), s));
    }
    for (int i = 0; i < AA.size(); i++)
    {
        Tcor.push_back(Treserve[i]);
    }
    return Tcor;
}

vector<double> TA2correction(TString name, int iter, vector<double> TT, vector<double> AA1, vector<double> AA2, int rbt = 10, int rbA = 12, double initialtL = -2e3, double initialtR = 3e3, double initialAL = -300, double initialAR = 300)
{
    bint = (initialtR - initialtL) / 0.5;
    cout << "bint: " << bint << endl;
    cout << "tL, tR: " << initialtL << "\t" << initialtR << endl;
    ht = new TH1D("ht", ";Time (ns);Counts", bint, initialtL, initialtR);
    htcor = new TH1D("htcor", ";Time (ns);Counts", bint, initialtL, initialtR);

    TH1D *hA[2];
    TH2D *hAT[2];
    hA[0] = new TH1D("hA0", ";A ;Counts", 3e3, initialAL, initialAR);
    hA[1] = new TH1D("hA1", ";A ;Counts", 3e3, initialAL, initialAR);
    hAT[0] = new TH2D("hAT0", "; A;TR (ns)", 3e3, initialAL, initialAR, bint, initialtL, initialtR);
    hAT[1] = new TH2D("hAT1", "; A;TR (ns)", 3e3, initialAL, initialAR, bint, initialtL, initialtR);

    TCanvas *cc;
    int CNum = 0;
    TH1F *htfit;
    TH1F *hAfit;
    TF1 *fitT;
    TF1 *fitA;
    TF1 *fitAT = new TF1("fitAT", "0", initialAL, initialAR);

    double Treserve[200000];
    double TAcor = 0;
    double Amean = 0;
    double Asigma = 0;
    double Tmean = 0;
    double Tsigma = 0;
    double fTL = 0;
    double fTR = 0;
    vector<double> Tcor;
    vector<double> AA[2];
    AA[0] = AA1;
    AA[1] = AA2;

    cc = cdC(CNum++);
    for (int s = 0; s < iter; s++)
    {
        for (int h = 0; h < 2; h++)
        {

            if (s != 0 || h != 0)
            {
                //if (Amean - 3 * Asigma <= 0)
                //    fitAT = profilefit(hAT, rbA, rbt * 4, Tmean - 6 * Tsigma, Tmean + 6 * Tsigma, 0.1, Amean + 3 * Asigma, Form("%s/%d", path.Data(), s-1));
                //else
                fitAT = profilefit(hAT[(h + 2 - 1) % 2], rbA * 2, rbt * 4, Tmean - 6 * Tsigma, Tmean + 6 * Tsigma, Amean - 3 * Asigma, Amean + 3 * Asigma, Form("%s_ch%dprofile%d", name.Data(), (h + 2 - 1) % 2, s - 1));
                if (!fitAT)
                {
                    cout << " the profilefit is failed! " << endl;
                    return Tcor;
                }
                ht->Reset();
                hAT[(h + 2 - 1) % 2]->Reset();
            }
            for (int i = 0; i < TT.size(); i++)
            {
                if (s == 0)
                    hA[h]->Fill(AA[h][i]);
                if (s == 0 && h == 0)
                {
                    //cout << "AA & TT: " << AA[i] << "\t" << TT[i] << endl;

                    Treserve[i] = TT[i];
                    fTL = initialtL;
                    fTR = initialtR;
                }
                else
                {

                    TAcor = fitAT->Eval(AA[(h + 2 - 1) % 2][i]);
                    //cout << "TAcor: " << TAcor << endl;
                    Treserve[i] = Treserve[i] - TAcor;
                    fTL = -1.5e3;
                    fTR = 1.5e3;
                }
                ht->Fill(Treserve[i]);
                hAT[h]->Fill(AA[h][i], Treserve[i]);
            }

            // * draw Tracklength
            cc = cdC(100);
            if (s == 0)
            {
                //DrawMyHist(hL, "", "", 1, 3);
                hA[h]->Draw();
                fitA = gausfit(hA[h], 100, 3, 3, rbA * 4, initialAL, initialAR);
                //return;
                if (!fitA)
                {
                    cout << "the hAfit hist is NULL " << endl;
                    return Tcor;
                }
                //fitA = (TF1 *)hAfit->GetFunction("fitU");
                Amean = fitA->GetParameter(1);
                Asigma = fitA->GetParameter(2);
                cout << "Amean=" << Amean << ",\tAsigma=" << Asigma << endl;
                cc->SaveAs(Form("%sch%dtot.png", name.Data(), h));
            }

            // * draw T
            cc = cdC(101);
            cc->Clear();
            //DrawMyHist(ht, "", "", 1, 3);
            ht->Draw();
            fitT = gausfit(ht, 100, 3, 3, rbt * 4, fTL, fTR);
            //fitT = gausfit(ht, 0.1, 3, 3, rbt * 20, fTL, fTR);
            if (!fitT)
            {
                cout << "the htfit hist is NULL " << endl;
                return Tcor;
            }
            //fitT = (TF1 *)htfit->GetFunction("fitU");
            Tmean = fitT->GetParameter(1);
            Tsigma = fitT->GetParameter(2);

            ht->SetNdivisions(505);
            sprintf(buff, "#sigma=%.0fps", Tsigma);
            TLatex *l = DrawMyLatex(buff, 0.2, 0.5);
            l->Draw();
            cout << "Tmean=" << Tmean << ",\tTsigma=" << Tsigma << endl;
            cc->SaveAs(Form("%s_ch%dTR_cor%d.png", name.Data(), h, s));

            /*
            //* draw T vs track
            cc = cdC(102);
            cc->Clear();
            //DrawMy2dHist(hAT, "", "", 1, 2);
            //ht->Rebin(2);
            hAT[h]->Draw("colz");
            cc->Modified();
            cc->Update();

            //fit = gausfit(ht, 20e-3, 3, 3, 1, fTL, fTR);
            //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
            //la = DrawMyLatex(buff, 0.2, 0.4);
            cc->SaveAs(Form("%s_ch%dhAT%d.png", name.Data(), h, s));
*/
        }
    }
    for (int i = 0; i < TT.size(); i++)
    {
        Tcor.push_back(Treserve[i]);
    }
    return Tcor;
}

void MatchTracker(TString fName1, TString fName2, TString fMatName, int force)
{

    //* check matched rootfile
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fMatName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fMatName << " is exist!" << endl;
            return;
        }
    }
    TFile *comFile = new TFile(fMatName, "recreate");
    if (!comFile->IsOpen())
    {
        cout << "ERROR! " << fMatName << " cant open!" << endl;
        return;
    }
    DefineRBCRhist();
    FTOFHit fFTOF;
    RBCRdata fCR;
    TTree *tree = new TTree("tree", "After matched");
    tree->Branch("FTOFdata", &fFTOF);
    tree->Branch("RBCRdata", &fCR);

    //* read FTOF data file
    FTOFHit *FTOFhit = new FTOFHit;
    TFile *FTOFfile = new TFile(fName1, "read");
    TTree *FTOFtree = (TTree *)FTOFfile->Get("tree");
    FTOFtree->SetBranchAddress("FTOFhit", &FTOFhit);

    //* read Tracker data file
    Int_t event;
    Double_t x0[4];
    Double_t y0[4];
    Double_t z0[4];
    Double_t x1[4];
    Double_t y1[4];
    Double_t z1[4];
    Bool_t sig_x0[4];
    Bool_t sig_y0[4];
    Bool_t sig_x1[4];
    Bool_t sig_y1[4];

    // List of branches
    TBranch *b_event;  //!
    TBranch *b_x0;     //!
    TBranch *b_y0;     //!
    TBranch *b_z0;     //!
    TBranch *b_x1;     //!
    TBranch *b_y1;     //!
    TBranch *b_z1;     //!
    TBranch *b_sig_x0; //!
    TBranch *b_sig_y0; //!
    TBranch *b_sig_x1; //!
    TBranch *b_sig_y1; //!
    TFile *trackerfile = new TFile(fName2, "read");
    TTree *fTree = (TTree *)trackerfile->Get("tree");
    fTree->SetMakeClass(1);
    fTree->SetBranchAddress("event", &event, &b_event);
    fTree->SetBranchAddress("x0", x0, &b_x0);
    fTree->SetBranchAddress("y0", y0, &b_y0);
    fTree->SetBranchAddress("z0", z0, &b_z0);
    fTree->SetBranchAddress("x1", x1, &b_x1);
    fTree->SetBranchAddress("y1", y1, &b_y1);
    fTree->SetBranchAddress("z1", z1, &b_z1);
    fTree->SetBranchAddress("sig_x0", sig_x0, &b_sig_x0);
    fTree->SetBranchAddress("sig_y0", sig_y0, &b_sig_y0);
    fTree->SetBranchAddress("sig_x1", sig_x1, &b_sig_x1);
    fTree->SetBranchAddress("sig_y1", sig_y1, &b_sig_y1);

    //* primary filter tracker data
    int NTracker = fTree->GetEntriesFast();
    vector<int> trackereventvec;
    vector<vector<double>> xposvec;
    vector<vector<double>> yposvec;
    vector<vector<double>> zposvec;
    vector<double> xpos;
    vector<double> ypos;
    vector<double> zpos;
    for (int i = 0; i < NTracker; i++)
    {
        xpos.clear();
        ypos.clear();
        zpos.clear();
        fTree->GetEntry(i);
        for (int s = 0; s < 4; s++)
        {
            if (x1[s] > -75 && x1[s] < 75 && y1[s] > -75 && y1[s] < 75)
            {
                xpos.push_back(x1[s]);
                ypos.push_back(y1[s]);
                zpos.push_back(z1[s]);
            }
        }
        if (zpos.size() > 2)
        {
            trackereventvec.push_back(event);
            xposvec.push_back(xpos);
            yposvec.push_back(ypos);
            zposvec.push_back(zpos);
        }
    }

    //* primary filter tracker data
    int NFTOF = FTOFtree->GetEntriesFast();
    int NT0ctr = 0;
    int NFTOFctr = 0;
    int NTrackerctr = 0;
    for (int i = 0, j = 0; i < NFTOF; i++)
    {
        DrawProcessbar(i, j, NFTOF);
        NT0ctr = 0;
        NFTOFctr = 0;
        NTrackerctr = 0;
        FTOFtree->GetEntry(i);
        for (int s = 0; s < 96; s++)
        {
            if (FTOFhit->lowtot[s] > 0)
                NFTOFctr++;
        }
        for (int s = 100; s < 104; s++)
        {
            if (FTOFhit->lowtot[s] > 0)
                NT0ctr++;
        }
        if (NFTOFctr > 0 && NT0ctr > 0)
        //if (NFTOFctr > 0 )
        //if (NT0ctr > 0 )
        {
            int iEntry = 0;
            iEntry = IsMatched(FTOFhit->event, trackereventvec);
            if (iEntry)
            {
                cout << FTOFhit->event << "\t" << iEntry << "\t" << trackereventvec[iEntry] << endl;
                NTrackerctr++;
                if (RebuildCRAngle(xposvec[iEntry], yposvec[iEntry], zposvec[iEntry], fCR, 0.9))
                {

                    fFTOF = *FTOFhit;
                    tree->Fill();
                    htheta->Fill(fCR.Mutheta);
                    hphi->Fill(fCR.Muphi);
                    hT0map->Fill(fCR.T0detX, fCR.T0detY);
                    hFTOFmap->Fill(fCR.FTOFdetX, fCR.FTOFdetY);
                }
            }
        }

        //* ToDo
    }
    comFile->WriteTObject(tree);
    TString path = GetPath(fMatName);
    DrawRBCRresults(Form("%s/RBCR.png", path.Data()));
    comFile->Flush();
    comFile->Close();
}
void RebuildT0(TString input, int force = 0)
{
    gStyle->SetOptFit(111);

    vector<double> TT, AA, tot1, tot2;
    TT.reserve(20e4);
    AA.reserve(20e4);
    bool flag = 0;
    TString filepath;
    TString rootname;
    rootname = GetFilename(input);
    filepath = GetPath(input);
    TString output = filepath + "/T0results.root";
    TFile *savefile = new TFile(output.Data(), "recreate");
    if (!savefile->IsOpen())
    {
        cout << "ERROR! " << output << " cant open!" << endl;
        return;
    }

    sprintf(buff, "%s/%s.dat", filepath.Data(), rootname.Data());
    if (!gSystem->AccessPathName(buff) && !force)
    {
        ifstream in;
        in.open(buff);
        double temp1 = 0, temp2 = 0;
        while (in && (!in.eof()))
        {

            in >> temp1 >> temp2;
            TT.push_back(temp1);
            AA.push_back(temp2);
        }
        flag = 1;
    }
    if (flag)
    {
        TAcorrection(filepath, 6, TT, AA);
        return;
    }
    ofstream out(buff);
    cout << "===> Create your data file: " << buff << endl;

    double Timestamp = 0;
    double Timestampcor = 0;
    double PMTtime = 0;
    double T0time[4] = {0};
    double T0timecor[4] = {0};
    double T0reTrack[4] = {0};
    double T0lowtot[4] = {0};
    double PMTtimecor = 0;
    double reTrack = 0;
    double reTrackSum = 0;
    int PMTcounter = 0;
    int validcnt = 0;
    double meanreTrack = 0;
    TTree *tree = new TTree("tree", "T0 results");
    tree->Branch("T0time", T0time, "T0time[4]/D");
    tree->Branch("T0lowtot", T0lowtot, "T0lowtot[4]/D");
    tree->Branch("T0timecor", T0timecor, "T0timecor[4]/D");
    tree->Branch("T0reTrack", T0reTrack, "T0reTrack[4]/D");

    cout << "Read rootfile: " << input.Data() << endl;
    TFile *f = new TFile(input.Data(), "READ");
    TTree *t = (TTree *)f->Get("tree");
    FTOFHit *fFTOF = new FTOFHit;
    RBCRdata *fCR = new RBCRdata;
    t->SetBranchAddress("FTOFdata", &fFTOF);
    t->SetBranchAddress("RBCRdata", &fCR);
    int N = t->GetEntriesFast();
    cout << "Total Matched events is: " << N << endl;
    for (int i = 0; i < N; i++)
    {
        t->GetEntry(i);
        //cout<<"> check <"<<endl;
        reTrackSum = 0;
        PMTtime = 0;

        PMTtimecor = 0;
        meanreTrack = 0;
        memset(T0time, 0, sizeof(T0time));
        memset(T0timecor, 0, sizeof(T0timecor));
        memset(T0reTrack, 0, sizeof(T0reTrack));
        PMTcounter = 0;
        for (int j = 100; j < 104; j++)
        {
            if (fFTOF->lowtot[j] > 0 && fFTOF->lowtot[j] < 2e3)
            {
                TVector3 Inpos(fCR->FTOFdetX, fCR->FTOFdetY, Zpos[T0id]);
                TVector3 Mudir(fCR->Mupx, fCR->Mupx, fCR->Mupz);
                if (j == 103)
                    reTrack = RebuildT0Track(Np, Inpos, Mudir, 3);
                else
                    reTrack = RebuildT0Track(Np, Inpos, Mudir, j - 100 + 1);
                //cout<<"reTrack: "<<reTrack<<endl;
                reTrackSum += reTrack;
                T0reTrack[j - 100] = reTrack;
                T0time[j - 100] = fFTOF->lowthtime[j];
                PMTtime += fFTOF->lowthtime[j];
                PMTtimecor += fFTOF->lowthtime[j] - reTrack / 298 * 1.5;
                T0timecor[j - 100] = fFTOF->lowthtime[j] - reTrack / 298 * 1.5;
                PMTcounter++;
                T0lowtot[j - 100] = fFTOF->lowtot[j];
            }
        }
        tree->Fill();
        if (T0time[0] && T0time[3])
        {
            Timestamp = T0time[0] - T0time[3];
            Timestampcor = T0timecor[0] - T0timecor[3];
            meanreTrack = T0reTrack[0] - T0reTrack[3];
            if (abs(Timestamp) < 1e3)
            {
                //if(1){
                TT.push_back(Timestamp);
                AA.push_back(meanreTrack);
                tot1.push_back(T0lowtot[0]);
                tot2.push_back(T0lowtot[3]);
                //cout << TT.back() << "\t" << AA.back() << endl;
                out << TT.back() << "\t" << AA.back() << endl;
            }
        }
    }
    savefile->WriteTObject(tree);

    TString name;
    vector<double> Tcor;

    name = filepath + "/TimeTrack";
    Tcor = TAcorrection(name, 6, TT, AA, 12, 10, -2e3, 3e3, -300, 300);

    name = filepath + "/TimeTOT";
    Tcor = TA2correction(name, 6, Tcor, tot1, tot2, 10, 10, -2e3, 2e3, 200, 1.5e3);

    /*
    name = filepath + "/TimeTOT1";
    Tcor = TAcorrection(name, 3, Tcor, tot1, 10, 10, -2e3, 2e3, 200, 1.5e3);
    name = filepath + "/TimeTOT2";
    TAcorrection(name, 3, Tcor, tot2, 10, 10, -2e3, 2e3, 200, 1.5e3);
*/
    savefile->Flush();
    //savefile->Close();
    //f->Close();
}
void Drawrisetime(TString fAnaName, int nbins = 100, double RL = 0, double RR = 2e3)
{
    cout << "open your file: " << fAnaName << endl;
    //TFile *hitFile = new TFile(fHitName, "recreate");
    TH1D *hrise[4];
    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[4];
    AGETHit *dethit;
    dethit = new AGETHit;
    TString fpath = GetPath(fAnaName);
    for (int s = 0; s < 4; s++)
    {
        cout << "Entry Detector " << s << " ..." << endl;
        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        tree[s]->SetBranchAddress("dethit", &dethit);
        hrise[s] = new TH1D(Form("hrise%d", s), Form("hrise%d", s), 200, RL, RR);
        int nentries = tree[s]->GetEntriesFast();
        cout << "Entries is " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {
            Long64_t ii = tree[s]->LoadTree(i);
            if (ii < 0)
                break;

            tree[s]->GetEntry(ii);
            if (dethit->wave_max > 7 * dethit->ped_rms && dethit->q > 5)
            {
                hrise[s]->Fill(dethit->risetime);
            }
        }
    }

    TCanvas *cc = cdC(0, 1400, 600);
    cc->Divide(4, 1);
    for (int s = 0; s < NumMM; s++)
    {
        cc->cd(s + 1);
        hrise[s]->Draw();
    }
    TString preName = fAnaName.Copy();
    preName = preName.Remove(preName.Length() - 5, 5);

    TString PngName;
    PngName = preName + "Trackerrisetime.png";
    cc->SaveAs(PngName);
}
void DrawDecodeInfo(TString fDecName, int nbins = 100)
{
    TH1F *hx[4];
    TH1F *hy[4];
    TH2F *hxy[4];
    TH1F *hampx[4];
    TH1F *hampy[4];
    TH1I *hstripx[4];
    TH1I *hstripy[4];
    int event;
    PositionData detector_data[NumMM];
    TBranch *b_dector[NumMM]; //!
    TFile *anaFile = new TFile(fDecName, "read");

    TTree *dtree = (TTree *)anaFile->Get("tree");
    dtree->SetMakeClass(1);
    dtree->SetBranchAddress("event", &event);
    for (int s = 0; s < NumMM; s++)
    {
        detector_data[s] = PositionData();
        string dec_name = "dector" + to_string(s);
        dtree->SetBranchAddress(dec_name.c_str(), &detector_data[s].sig, &b_dector[s]);

        hx[s] = new TH1F(Form("hx%d", s), "hit position in x axis", nbins, 0, 400);
        hy[s] = new TH1F(Form("hy%d", s), "hit position in y axis", nbins, 0, 400);
        hxy[s] = new TH2F(Form("hxy%d", s), "hit position", nbins / 2, 0, 400, nbins / 2, 0, 400);
        hampx[s] = new TH1F(Form("hampx%d", s), "total charge of hit in x axis", nbins, 0, 4e3);
        hampy[s] = new TH1F(Form("hampy%d", s), "total charge of hit in y axis", nbins, 0, 4e3);
        hstripx[s] = new TH1I(Form("hstripx%d", s), "number of fired asic channels in x axis", 20, 0, 20);
        hstripy[s] = new TH1I(Form("hstripy%d", s), "number of fired asic channels in y axis", 20, 0, 20);
    }
    int nentries = dtree->GetEntriesFast();
    cout << " --> nentries= " << nentries << endl;
    for (int i = 0; i < nentries; i++)
    {

        dtree->GetEntry(i);
        for (int s = 0; s < NumMM; s++)
        {
            if (detector_data[s].hit_strip_num_x > 0)
            {

                hx[s]->Fill(detector_data[s].x);
                hampx[s]->Fill(detector_data[s].hit_amp_x);
                hstripx[s]->Fill(detector_data[s].hit_strip_num_x);
            }

            if (detector_data[s].hit_strip_num_y > 0)
            {
                hy[s]->Fill(detector_data[s].y);
                hampy[s]->Fill(detector_data[s].hit_amp_y);
                hstripy[s]->Fill(detector_data[s].hit_strip_num_y);
            }
            if (detector_data[s].hit_strip_num_y > 0 && detector_data[s].hit_strip_num_x > 0)
                hxy[s]->Fill(detector_data[s].x, detector_data[s].y);
        }
    }
    TCanvas *cc = cdC(0, 1400, 1000);
    cc->Divide(4, 3);
    for (int s = 0; s < NumMM; s++)
    {
        cc->cd(s + 1);
        hx[s]->Draw();

        cc->cd(s + 5);
        hy[s]->Draw();
        cc->cd(s + 9);
        hxy[s]->Draw("colz");
    }
    TString preName = fDecName.Copy();
    preName = preName.Remove(preName.Length() - 5, 5);

    TString PngName;
    PngName = preName + "TrackerdecMap.png";
    cc->SaveAs(PngName);

    cc = cdC(1, 1600, 1000);
    cc->Divide(4, 4);
    TLatex *la;
    for (int s = 0; s < NumMM; s++)
    {
        cc->cd(s + 1);
        hampx[s]->Draw();
        hampx[s]->Fit("landau");
        double peakx = hampx[s]->GetFunction("landau")->GetParameter(1);
        la = DrawMyLatex(Form("MPV=%.f", peakx), 0.64, 0.5);
        la->Draw("same");

        cc->cd(s + 5);
        hampy[s]->Draw();
        hampy[s]->Fit("landau");
        double peaky = hampy[s]->GetFunction("landau")->GetParameter(1);
        la = DrawMyLatex(Form("MPV=%.f", peaky), 0.64, 0.5);
        la->Draw("same");
        cc->cd(s + 9);
        hstripx[s]->Draw();
        cc->cd(s + 13);
        hstripy[s]->Draw();
    }
    preName = fDecName.Copy();
    preName = preName.Remove(preName.Length() - 5, 5);

    PngName = preName + "TrackerdecTotalamp.png";
    cc->SaveAs(PngName);
}
#if 1
void PrintTracker(TString fAnaName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fAnaName, fStat);
        if (fStat.fSize == 0)
        {
            cout << "ERROR! " << fAnaName << " isn't exist!" << endl;
            return;
        }
    }

    int chNo = 64; // 64*6=384
    //TFile *hitFile = new TFile(fHitName, "recreate");
    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree[4];
    AGETHit *dethit;
    dethit = new AGETHit;
    //TBranch *b_dethit;

    TH1I *hDeteff[4];
    TH2I *hDetHit[4];
    TH1I *hXid[4];
    TH1I *hYid[4];
    TH1I *hXHit[4];
    TH1I *hYHit[4];
    vector<AGETHit> xhit[4];
    vector<AGETHit> yhit[4];

    TH1I *hid = new TH1I("hid", "hid", 100e4, 1, 100e4 + 1);
    for (int s = 0; s < 4; s++)
    {

        tree[s] = (TTree *)anaFile->Get(Form("tree%d", s));
        tree[s]->SetBranchAddress("dethit", &dethit);

        hXid[s] = new TH1I(Form("hXid%d", s), Form("hXid%d", s), 100e4, 1, 100e4 + 1);
        hYid[s] = new TH1I(Form("hYid%d", s), Form("hYid%d", s), 100e4, 1, 100e4 + 1);
        //x-2,y-1
        hDeteff[s] = new TH1I(Form("hDeteff%d", s), "", 3, 0, 3);
        //leg->AddEntry(hDeteff[i], Form("Det%d",i), "lp");

        hDetHit[s] = new TH2I(Form("hDetHit%d", s), "", chNo, 0, chNo, chNo, 0, chNo);
        hXHit[s] = new TH1I(Form("HitX%d", s), "X Strip hits for single event", chNo, 0, chNo);
        hYHit[s] = new TH1I(Form("HitY%d", s), "Y Strip hits for single event", chNo, 0, chNo);
        //memset(xhit, 0, sizeof(xhit));
        //memset(yhit, 0, sizeof(yhit));

        int nentries = tree[s]->GetEntriesFast();
        cout << " --> nentries= " << nentries << endl;
        for (int i = 0; i < nentries; i++)
        {

            Long64_t ii = tree[s]->LoadTree(i);
            if (ii < 0)
                break;
            if (ii % 10000 == 0)
                cout << " --> Entry$= " << ii << endl;
            tree[s]->GetEntry(ii);
            hid->Fill(dethit->event);
            //if (dethit->wave_max > 7 * dethit->ped_rms && dethit->q > 5)
            if (1)
            {
                // x axis has a hit
                if (dethit->chn[0] != -999 && dethit->wave_max > 7 * dethit->ped_rms)
                {

                    hXHit[s]->Fill(dethit->chn[0]);
                    xhit[s].push_back(*dethit);
                    hXid[s]->Fill(dethit->event);
                }

                // y axis has a hit
                if (dethit->chn[1] != -999 && dethit->wave_max > 6 * dethit->ped_rms)
                {

                    hYHit[s]->Fill(dethit->chn[1]);
                    yhit[s].push_back(*dethit);
                    hYid[s]->Fill(dethit->event);
                }
            }
        }
        cout << "MM" << s << " read completely" << endl;
    }
    //hid->Draw();
    //return;
    RecordStart = GetFirstTrig(hid);
    RecordEnd = GetLastTrig(hid);
    hid->GetXaxis()->SetRangeUser(RecordStart - 5, RecordEnd + 5);
    vector<int> xtemp;
    vector<int> ytemp;

    Color_t clr[] = {2, 4, 1};
    Style_t stl[] = {7, 7, 1};
    TLegend *leg = DrawMyLeg(0.7, 0.4, 0.85, 0.55);
    TLatex *la;
    TCanvas *c;

    for (int s = RecordStart, j = 0; s < RecordEnd + 1; s++)
    {
        DrawProcessbar(s - RecordStart, j, RecordEnd + 1 - RecordStart);
        for (int i = 0; i < 4; i++)
        {
            xtemp.clear();
            ytemp.clear();
            int nsize = xhit[i].size();
            for (int j = 0; j < nsize; j++)
            {
                if (xhit[i].at(j).event == s)
                    xtemp.push_back(xhit[i].at(j).chn[0]);
            }
            nsize = yhit[i].size();
            for (int j = 0; j < nsize; j++)
            {

                if (yhit[i].at(j).event == s)
                    ytemp.push_back(yhit[i].at(j).chn[1]);
            }

            // Fill 2d-hit map
            /*
            for (int k = 0; k < (int)xtemp.size(); k++)
                for (int j = 0; j < (int)ytemp.size(); j++)
                {
                    hDetHit[i]->Fill(xtemp.at(k), ytemp.at(j));
                }
            */
            if (xtemp.size() > 0)
            {
                hDeteff[i]->Fill(2);
            }
            if (ytemp.size() > 0)
            {

                hDeteff[i]->Fill(1);
            }
        }
    }
    cout << "Hist fill completely" << endl;

    c = new TCanvas("c3", "c3", 1500, 1200);
    c->Divide(4, 4);
    double detEff[2];
    TString XY[] = {"Y", "X"};
    for (int i = 0; i < 4; i++)
    {
        c->cd(i + 1);
        detEff[0] = 0;
        detEff[1] = 0;
        hDeteff[i]->Draw();
        for (int j = 0; j < 2; j++)
        {

            detEff[j] = hDeteff[i]->GetBinContent(j + 2) / NumofTriggers;
            sprintf(buff, "MM%d%s,Eff=%.1f%%", i, XY[j].Data(), detEff[j] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * j, 42, 0.05);
            la->Draw("same");
        }
        c->cd(5 + i);
        hXHit[i]->Draw();
        c->cd(9 + i);
        hYHit[i]->Draw();
        c->cd(13 + i);
        //hDetHit[i]->SetMarkerStyle(4);
        hDetHit[i]->Draw("colz");
    }
    TString PngName;
    PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "HitDistributed.png";
    c->SaveAs(PngName);
    //hitFile->Close();
    //anaFile->Close();
    c = new TCanvas("c4", "c4", 1500, 1200);
    c->Divide(2, 4);
    for (int i = 0; i < 4; i++)
    {
        c->cd(i * 2 + 1);
        hXid[i]->Draw();
        hXid[i]->GetXaxis()->SetRangeUser(GetFirstTrig(hXid[i]) - 5, GetLastTrig(hXid[i]) + 5);
        c->cd(i * 2 + 2);
        hYid[i]->Draw();
        hYid[i]->GetXaxis()->SetRangeUser(GetFirstTrig(hYid[i]) - 5, GetLastTrig(hYid[i]) + 5);
    }
    PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "EventIDDistributed.png";
    c->SaveAs(PngName);

    //anaFile->Close();
    cout << "--> Efficiency results has been calculated ! " << endl;
}
#endif
void PrintT0(TString fileName, TString opfileName, int force = 0)
{

    int N = sizeof(T0ELEchid) / sizeof(T0ELEchid[0]);

    TFile *opfile;
    FileStat_t fStat;
    gSystem->GetPathInfo(opfileName, fStat);
    if (fStat.fSize != 0 && force != 1)
    {

        cout << "ERROR! " << opfileName << " is exist!" << endl;
        opfile = new TFile(opfileName, "read");
    }
    else
        opfile = new TFile(opfileName, "recreate");
    cout << "====>>  Create the root file : " << opfileName << endl;

    TGaxis::SetMaxDigits(3);

    Color_t clr[] = {8, 2, 4, 1};                                     //color
    Style_t fst[] = {1001, 3001, 3002, 3003, 3004, 3005, 3006, 3007}; //

    //sprintf(path, "/mnt/f/CRsys/CosmicRayTest/%s%s", date, root);
    //sprintf(path, "/mnt/d/ExpDATA/labtest/CosmicRay/%s%s", date, root);

    //sprintf(rootname,"%s/beamtest%s/Analysismatch%s",rootpath,date,date);

    //double LowThTimediff;
    double tL = -500e3, tR = -300e3;
    double rtL = -1e3, rtR = 0.2e3;
    double TcutL = 1;
    double TcutR = 2e3; //reffered Tube;
    //TcutL=1;
    //TcutR=4e3; //reffered Tube;

    double UL = 0, UR = 30e3;
    double UcutL = 1, UcutR = 3e3; //reffered Tube;

    //UcutL=1;
    //UcutR=4e3; //reffered Tube;

    //int rbt=10;
    //int rbU=8;
    TH1D *h1d;
    TH1D *h1dpick;
    TH1D *h2d;
    TH1D *h2dpick;
    TFile *file = new TFile(fileName, "read");
    //TTree* tree =  (TTree*)file->Get("Analysis");
    TTree *tree = (TTree *)file->Get("tree");
    cout << "====>>  Start to open the root file : " << fileName << endl;

    FTOFHit *FTOFhit = new FTOFHit;
    tree->SetBranchAddress("FTOFhit", &FTOFhit);

    TH1D *hLowt = new TH1D("hLowt", "Time Resolution;LowTh time (ps);Counts", 500, tL, tR);
    TH1D *hLowt_Pick = (TH1D *)hLowt->Clone("hLowt_Pick");

    TH1D *hHight = new TH1D("hHight", "Time Resolution;HighTh time (ps);Counts", 500, tL, tR);
    TH1D *hHight_Pick = (TH1D *)hHight->Clone("hHight_Pick");

    TH1D *hRise = new TH1D("hRise", "Rise Time; Rise time (ps);Counts", 0.2e3, rtL, rtR);
    TH1D *hRise_Pick = (TH1D *)hRise->Clone("hRise_Pick");

    TH1D *hHighU = new TH1D("hHighU", "HighTOT;HighTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hLowU = new TH1D("hLowU", "LowTOT;LowTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hHighU_Pick = (TH1D *)hHighU->Clone("hHighU_Pick");
    TH1D *hLowU_Pick = (TH1D *)hLowU->Clone("hLowU_Pick");
    TH2D *h2dU = new TH2D("h2dU", ";LowTOT (ps)T;HighTOT (ps)", 500, UL, UR, 500, UL, UR);

    TH1D *hPMTID = new TH1D("hPMTID", ";PMTID;Counts", 4, 0, 4);
    TH1D *hPMTID_Pick = (TH1D *)hPMTID->Clone("hPMTID_Pick");

    TH1D *hNPMT = new TH1D("hNPMT", ";Number of Fired PMT;Counts", 4 + 1, 0, 4 + 1);
    TH1D *hNPMT_Pick = (TH1D *)hNPMT->Clone("hNPMT_Pick");

    TH2D *h2dPMT = new TH2D("h2dPMT", ";Number of Fired PMT;Timediff (ps)", 4 + 1, 0, 4 + 1, 100, tL, tR);
    TH2D *h2dPMT_Pick = (TH2D *)h2dPMT->Clone("h2dPMT_Pick");
    TH2D *hHighAT = new TH2D(buff, "hHighAT;HighTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowAT = new TH2D(buff, "hLowAT;LowTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowRT = new TH2D(buff, "hLowRT;LowTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);
    TH2D *hHighRT = new TH2D(buff, "hHighRT;HighTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);

    TH1D *ho[100];
    TH1D *hho[100];

    TCanvas *c;
    TCanvas *co;
    TLegend *leg;
    TLegend *lego = DrawMyLeg();
    TLatex *ll;

    int NEvents = tree->GetEntries();
    int IDIndex = -1;
    int eventcounter = 0;
    int counter = 0;

    vector<int> TOTIDIndex;
    TOTIDIndex.clear();
    for (int j = 0; j < N; j++)
    {

        bool LoopOn = 1;
        hLowt->Reset();
        hHight->Reset();
        hRise->Reset();
        hLowU->Reset();
        hHighU->Reset();
        hLowAT->Reset();
        hHighAT->Reset();
        hLowRT->Reset();
        hHighRT->Reset();
        h2dU->Reset();

        hLowt_Pick->Reset();
        hHight_Pick->Reset();
        hRise_Pick->Reset();
        hLowU_Pick->Reset();
        hHighU_Pick->Reset();

        tree->GetEntry(0);

        for (int i = 0; i < NEvents; i++)
        {
            //for(int i=0; i<100; i++){
            tree->GetEntry(i);
            if (j == 0)
            {
                double PMTtime = 0;
                double PMTtime_Pick = 0;
                int PMTcounter = 0;
                int PMTcounter_Pick = 0;
                for (int k = 0; k < N; k++)
                {
                    if (FTOFhit->lowtot[T0ELEchid[k]] > UcutL)
                    {
                        PMTtime += FTOFhit->triggertimediff[k];
                        PMTcounter++;
                        hPMTID->Fill(T0Detchid[k]);
                        if (FTOFhit->lowtot[T0ELEchid[k]] < UcutR)
                        {
                            PMTcounter_Pick++;
                            PMTtime_Pick += FTOFhit->triggertimediff[T0ELEchid[k]];
                            hPMTID_Pick->Fill(T0Detchid[k]);
                        }
                    }
                }
                if (PMTcounter > 0)
                {
                    eventcounter++;
                    hNPMT->Fill(PMTcounter);
                    h2dPMT->Fill(PMTcounter, PMTtime / PMTcounter);
                    hNPMT_Pick->Fill(PMTcounter_Pick);
                    h2dPMT_Pick->Fill(PMTcounter_Pick, PMTtime_Pick / PMTcounter_Pick);
                }
            }

            //else
            //LoopSwitch = 0;
            double LowThTimediff;
            double HighThTimediff;
            double risetime;
            if (FTOFhit->lowtot[T0ELEchid[j]] > UcutL)
            {
                LowThTimediff = FTOFhit->triggertimediff[T0ELEchid[j]];
                HighThTimediff = FTOFhit->highthtime[T0ELEchid[j]] + FTOFhit->triggertimediff[T0ELEchid[j]] - FTOFhit->lowthtime[T0ELEchid[j]];
                risetime = FTOFhit->highthtime[T0ELEchid[j]] - FTOFhit->lowthtime[T0ELEchid[j]];
                hLowt->Fill(LowThTimediff);
                hHight->Fill(HighThTimediff);
                hRise->Fill(risetime);
                hLowU->Fill(FTOFhit->lowtot[T0ELEchid[j]]);
                hHighU->Fill(FTOFhit->hightot[T0ELEchid[j]]);
                hLowAT->Fill(FTOFhit->lowtot[T0ELEchid[j]], LowThTimediff);
                hHighAT->Fill(FTOFhit->hightot[T0ELEchid[j]], LowThTimediff);
                hLowRT->Fill(FTOFhit->lowtot[T0ELEchid[j]], risetime);
                hHighRT->Fill(FTOFhit->hightot[T0ELEchid[j]], risetime);
                h2dU->Fill(FTOFhit->lowtot[T0ELEchid[j]], FTOFhit->hightot[T0ELEchid[j]]);

                //if (LowTOT[T0ELEchid[j]] < UcutR && LowTOT[T0ELEchid[j]] > UcutL && (LowTOT[T0ELEchid[j]] > LowTOT[other[0]] || LowTOT[T0ELEchid[j]] > LowTOT[other[1]] || LowTOT[T0ELEchid[j]] > LowTOT[other[2]]))
                if (FTOFhit->lowtot[T0ELEchid[j]] < UcutR)

                {
                    hLowt_Pick->Fill(LowThTimediff);
                    hHight_Pick->Fill(HighThTimediff);
                    hRise_Pick->Fill(risetime);
                    hLowU_Pick->Fill(FTOFhit->lowtot[T0ELEchid[j]]);
                    hHighU_Pick->Fill(FTOFhit->hightot[T0ELEchid[j]]);
                }
            }
        }
        ValidT0 = eventcounter;
        cout << "check !" << endl;
        TString PngName = fileName.Copy();
        PngName = PngName.Remove(PngName.Length() - 5, 5) + "T" + TString::Itoa(T0Detchid[j], 10);
        TString cName = Form("T%d", T0Detchid[j]);

        sprintf(buff, "T%d", T0Detchid[j]);
        ll = DrawMyLatex(buff);
        //
        // ---------draw LowTOT vs HighTOT --------//
        //
        c = cdC(counter++);
        //ht->Rebin(2);
        h2dU->Draw("colz");
        h2dU->GetXaxis()->SetNdivisions(505);
        DrawMy2dHist(h2dU, "", "", 1, 2);
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        ll->Draw();
        //sprintf(buff, "%s2dU.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_2dU", cName.Data()));
        //return;

        //
        // *** Low Threshold Time distribution ***
        // **
        c = cdC(counter++);
        hLowt->Draw();
        DrawMyHist(hLowt, "", "", 1, 2, 1);
        DrawMyHist(hLowt_Pick, "", "", 2, 2, 7);
        hLowt_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hLowt, "No cut", "l");
        leg->AddEntry(hLowt_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sLowThTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowThTime", cName.Data()));

        //
        // *** HighThreshold Time distribution ***
        // **
        c = cdC(counter++);
        //c[0]->SetLogy();
        hHight->Draw();
        DrawMyHist(hHight, "", "", 1, 2, 1);
        DrawMyHist(hHight_Pick, "", "", 2, 2, 7);
        hHight_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hHight, "No cut", "l");
        leg->AddEntry(hHight_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sHighThTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighThTime", cName.Data()));

        //
        // *** Rise Time distribution ***
        // *********
        c = cdC(counter++);
        hRise->Draw();
        DrawMyHist(hRise, "", "", 1, 2, 1);
        DrawMyHist(hRise_Pick, "", "", 2, 2, 7);
        hRise_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hRise, "No cut", "l");
        leg->AddEntry(hRise_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();

        //sprintf(buff, "%sRiseTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_RiseTime", cName.Data()));

        //
        // *** Low TOT distribution ***
        // *********
        c = cdC(counter++);
        hLowU->GetXaxis()->SetRangeUser(UL, UR);
        hLowU->Draw();

        DrawMyHist(hLowU, "", "", 1, 2, 1);
        DrawMyHist(hLowU_Pick, "", "", 2, 2, 7);
        hLowU_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hLowU, "No cut", "l");
        leg->AddEntry(hLowU_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sLowTOTL.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowTOTL", cName.Data()));

        hLowU->GetXaxis()->SetRangeUser(0, 1.5e3);
        //sprintf(buff, "%sLowTOTS.png", PngName.Data());
        //c->SaveAs(buff);
        //opfile->WriteObject(c, Form("%s_LowTOTS", cName.Data()));

        sprintf(buff, "ho%d", T0Detchid[j]);
        TOTIDIndex.push_back(T0Detchid[j]);
        ho[T0Detchid[j]] = (TH1D *)hLowU->Clone(buff);

        //
        // *** High TOT distribution ***
        // *********
        c = cdC(counter++);
        hHighU->GetXaxis()->SetRangeUser(UL, UR);
        hHighU->Draw();

        DrawMyHist(hHighU, "", "", 1, 2, 1);
        DrawMyHist(hHighU_Pick, "", "", 2, 2, 7);
        hHighU_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hHighU, "No cut", "l");
        leg->AddEntry(hHighU_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sHighTOTL.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighTOTL", cName.Data()));

        hHighU->GetXaxis()->SetRangeUser(0, 1.5e3);
        //sprintf(buff, "%sHighTOTS.png", PngName.Data());
        //c->SaveAs(buff);

        sprintf(buff, "hho%d", T0Detchid[j]);
        hho[T0Detchid[j]] = (TH1D *)hHighU->Clone(buff);

        //
        // *** LowTOT distribution - Risetime ***
        // *********
        c = cdC(counter++);
        hLowRT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hLowRT, "", "");

        ll->Draw();
        //sprintf(buff, "%sLowRT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowRT", cName.Data()));

        //
        // *** HighTOT distribution - Risetime ***
        // *********
        c = cdC(counter++);
        hHighRT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hHighRT, "", "");

        ll->Draw();
        //sprintf(buff, "%sHighRT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighRT", cName.Data()));

        //
        // *** LowThTime-LowTOT distribution ***
        // *********
        c = cdC(counter++);
        hLowAT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hLowAT, "", "");

        ll->Draw();
        //sprintf(buff, "%sLowAT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowAT", cName.Data()));

        //
        // *** LowThTime-HighTOT distribution ***
        // *********
        c = cdC(counter++);
        hHighAT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);

        DrawMy2dHist(hHighAT, "", "");

        ll->Draw();
        //sprintf(buff, "%sHighAT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighAT", cName.Data()));

        cout << "loop state:" << LoopOn << ", No. " << j << endl;
        if (!LoopOn)
            break;
    }

    TString PngName = fileName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5);

#if 0

    //
    // ---------draw tot of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 3);
    TCanvas *p1;
    TCanvas *p2;
    TCanvas *p3;
    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowTOTL", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighTOTL", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_2dU", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
    }
    sprintf(buff, "%sTOT_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw time of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 3);

    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowThTime", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighThTime", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_RiseTime", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
    }
    sprintf(buff, "%sTime_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw 2d relationship  of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 4);

    TCanvas *p4;
    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowRT", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighRT", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_LowAT", T0Detchid[i]));
        p4 = (TCanvas *)opfile->Get(Form("T%d_HighAT", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 1 * N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 3 * N);
        p4->DrawClonePad();
    }
    sprintf(buff, "%s2Drelationship_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw NPMT vs TR --------//
    //
    c = cdC(counter++);
    DrawMy2dHist(h2dPMT, "", "", 1, 2);
    //ht->Rebin(2);
    h2dPMT->Draw("colz");
    h2dPMT->GetXaxis()->SetNdivisions(505);
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    //sprintf(buff, "%s2dPMT.png", PngName.Data());
    //c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_T0_vs_Time"));
    //return;

    

    //

#endif
    //
    // ---------draw efficiency of PMT--------//
    //
    double PMTEff[N];
    double PMTEff_Pick[N];
    c = cdC(counter++);
    DrawMyHist(hPMTID, "PMT ID", "Counts", 1, 3);
    //ht->Rebin(2);
    hPMTID->Draw();
    hPMTID->GetXaxis()->SetNdivisions(505);
    c->Update();
    SetEffstats(c, hPMTID, NumofTriggers, RecordTriggers, ValidT0);
    //DrawMyHist(hPMTID_Pick, "", "", 2, 2, 7);
    //hPMTID_Pick->Draw("same");
    //leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    //leg->AddEntry(hPMTID, "No cut", "l");
    //leg->AddEntry(hPMTID_Pick, "With cut", "l");
    //leg->Draw();
    for (int i = 0; i < N; i++)
    {

        PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / RecordTriggers;
        //PMTEff_Pick[i] = hPMTID_Pick->GetBinContent(hPMTID_Pick->FindBin(i)) / RecordTriggers;
        sprintf(buff, "PMTID=%d,Eff=%.1f%%", i, PMTEff[i] * 100);
        ll = DrawMyLatex(buff, 0.3, 0.7 - 0.07 * i, 42, 0.05, 2);
        ll->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_PMT.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_PMT"));

    //return;

    //
    // *** Number of fired PMT ***
    // **
    double Eff[N + 1];
    double Eff_Pick[N + 1];
    TLatex *la;

    c = cdC(counter++);

    //DrawMyHist(hNPMT_Pick, "", "", 2, 2, 7);
    //hNPMT_Pick->Draw("SAME");
    //leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    //leg->AddEntry(hNPMT, "No cut", "l");
    //leg->AddEntry(hNPMT_Pick, "With cut", "l");
    //leg->Draw();
    int noevent = hNPMT->GetBinContent(hNPMT->FindBin(0));
    if (noevent < 10)
    {
        noevent = RecordTriggers - ValidT0;
        hNPMT->SetBinContent(hNPMT->FindBin(0), noevent);
    }

    DrawMyHist(hNPMT, "Num of fired PMT", "Counts", 1, 3);
    //ht->Rebin(2);
    hNPMT->Draw();
    hNPMT->GetXaxis()->SetNdivisions(505);
    c->Update();
    SetEffstats(c, hNPMT, NumofTriggers, RecordTriggers, ValidT0);
    Eff[0] = hNPMT->GetBinContent(hNPMT->FindBin(0)) / RecordTriggers;
    sprintf(buff, "No PMT response, Ratio=%.1f%%", Eff[0] * 100);
    la = DrawMyLatex(buff, 0.25, 0.7, 42, 0.05, 1);
    la->Draw("same");
    for (int i = 1; i < N + 1; i++)
    {
        Eff[i] = hNPMT->GetBinContent(hNPMT->FindBin(i)) / RecordTriggers;
        //Eff_Pick[i] = hNPMT_Pick->GetBinContent(hNPMT_Pick->FindBin(i)) / eventcounter;
        sprintf(buff, "NPMT=%d,Ratio=%.1f%%,#color[2]{(Valid)%.1f%%}", i, Eff[i] * 100, Eff[i] * 100 / ValidT0 * RecordTriggers);
        la = DrawMyLatex(buff, 0.25, 0.7 - 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_T0.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_T0"));

    //
    // ---------draw LowTOT together --------//
    //
    int temp;
    int max_y_content;
    sort(TOTIDIndex.begin(), TOTIDIndex.end());
    co = cdC(counter++);
    temp = 0;
    max_y_content = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        FillHistArea(ho[TOTIDIndex.at(i)], clr[i], 0.8, fst[2]);
        sprintf(buff, "T%d", TOTIDIndex.at(i));
        lego->AddEntry(ho[TOTIDIndex.at(i)], buff, "lpf");
        temp = ho[TOTIDIndex.at(i)]->GetMaximum();
        max_y_content = max_y_content > temp ? max_y_content : temp;
        if (i == 0)
            ho[TOTIDIndex.at(i)]->Draw();
        else
            ho[TOTIDIndex.at(i)]->Draw("same");
    }
    ho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max_y_content * 1.1);
    ho[TOTIDIndex.at(0)]->SetStats(0);
    co->Modified();
    lego->Draw();
    sprintf(buff, "%sLowTOTDrawtogether.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_lowtot"));

    //
    // ---------draw HighTOT together --------//
    //
    lego->Clear();
    co = cdC(counter++);
    temp = 0;
    max_y_content = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        FillHistArea(hho[TOTIDIndex.at(i)], clr[i], 0.1, fst[2]);
        sprintf(buff, "T%d", TOTIDIndex.at(i));
        lego->AddEntry(hho[TOTIDIndex.at(i)], buff, "lpf");
        temp = hho[TOTIDIndex.at(i)]->GetMaximum();
        max_y_content = max_y_content > temp ? max_y_content : temp;

        if (i == 0)
            hho[TOTIDIndex.at(i)]->Draw();
        else
            hho[TOTIDIndex.at(i)]->Draw("same");
    }
    hho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max_y_content * 1.1);
    hho[TOTIDIndex.at(0)]->SetStats(0);
    co->Modified();
    lego->Draw();
    sprintf(buff, "%sHighTOTDrawtogether.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_hightot"));
}
/*
void PrintRawT0(TString fileName, TString opfileName, int force = 0)
{

    int N = sizeof(T0ELEchid) / sizeof(T0ELEchid[0]);

    TFile *opfile;
    FileStat_t fStat;
    gSystem->GetPathInfo(opfileName, fStat);
    if (fStat.fSize != 0 && force != 1)
    {

        cout << "ERROR! " << opfileName << " is exist!" << endl;
        opfile = new TFile(opfileName, "read");
    }
    else
        opfile = new TFile(opfileName, "recreate");
    cout << "====>>  Create the root file : " << opfileName << endl;

    TGaxis::SetMaxDigits(3);

    Color_t clr[] = {8, 2, 4, 6};                                     //color
    Style_t fst[] = {1001, 3001, 3002, 3003, 3004, 3005, 3006, 3007}; //

    //sprintf(path, "/mnt/f/CRsys/CosmicRayTest/%s%s", date, root);
    //sprintf(path, "/mnt/d/ExpDATA/labtest/CosmicRay/%s%s", date, root);

    //sprintf(rootname,"%s/beamtest%s/Analysismatch%s",rootpath,date,date);

    //double LowThTimediff;
    double tL = -800e3, tR = -600e3;
    double rtL = -1e3, rtR = 0.2e3;
    double TcutL = 1;
    double TcutR = 2e3; //reffered Tube;
    //TcutL=1;
    //TcutR=4e3; //reffered Tube;

    double UL = 0, UR = 30e3;
    double UcutL = 1, UcutR = 3e3; //reffered Tube;

    //UcutL=1;
    //UcutR=4e3; //reffered Tube;

    //int rbt=10;
    //int rbU=8;
    TH1D *h1d;
    TH1D *h1dpick;
    TH1D *h2d;
    TH1D *h2dpick;
    TFile *file = new TFile(fileName, "read");
    //TTree* tree =  (TTree*)file->Get("Analysis");
    TTree *tree = (TTree *)file->Get("tree");
    cout << "====>>  Start to open the root file : " << fileName << endl;

    int TrackID;
    int ChID;
    double LowThTime;
    double HighThTime;
    double LowTOT;
    double HighTOT;
    tree->SetBranchAddress("TrackID", &TrackID);
    tree->SetBranchAddress("ChID", &ChID);
    tree->SetBranchAddress("LowThTime", &LowThTime);
    tree->SetBranchAddress("HighThTime", &HighThTime);
    tree->SetBranchAddress("LowTOT", &LowTOT);
    tree->SetBranchAddress("HighTOT", &HighTOT);

    TH1D *hLowt = new TH1D("hLowt", "Time Resolution;LowTh time (ps);Counts", 500, tL, tR);
    TH1D *hLowt_Pick = (TH1D *)hLowt->Clone("hLowt_Pick");

    TH1D *hHight = new TH1D("hHight", "Time Resolution;HighTh time (ps);Counts", 500, tL, tR);
    TH1D *hHight_Pick = (TH1D *)hHight->Clone("hHight_Pick");

    TH1D *hRise = new TH1D("hRise", "Rise Time; Rise time (ps);Counts", 0.2e3, rtL, rtR);
    TH1D *hRise_Pick = (TH1D *)hRise->Clone("hRise_Pick");

    TH1D *hHighU = new TH1D("hHighU", "HighTOT;HighTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hLowU = new TH1D("hLowU", "LowTOT;LowTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hHighU_Pick = (TH1D *)hHighU->Clone("hHighU_Pick");
    TH1D *hLowU_Pick = (TH1D *)hLowU->Clone("hLowU_Pick");
    TH2D *h2dU = new TH2D("h2dU", ";LowTOT (ps)T;HighTOT (ps)", 500, UL, UR, 500, UL, UR);

    TH1D *hPMTID = new TH1D("hPMTID", ";PMTID;Counts", 4, 0, 4);
    TH1D *hPMTID_Pick = (TH1D *)hPMTID->Clone("hPMTID_Pick");

    TH1D *hNPMT = new TH1D("hNPMT", ";Number of Fired PMT;Counts", 4 + 1, 0, 4 + 1);
    TH1D *hNPMT_Pick = (TH1D *)hNPMT->Clone("hNPMT_Pick");

    TH2D *h2dPMT = new TH2D("h2dPMT", ";Number of Fired PMT;Timediff (ps)", 4 + 1, 0, 4 + 1, 100, tL, tR);
    TH2D *h2dPMT_Pick = (TH2D *)h2dPMT->Clone("h2dPMT_Pick");
    TH2D *hHighAT = new TH2D(buff, "hHighAT;HighTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowAT = new TH2D(buff, "hLowAT;LowTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowRT = new TH2D(buff, "hLowRT;LowTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);
    TH2D *hHighRT = new TH2D(buff, "hHighRT;HighTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);

    TH1D *ho[100];
    TH1D *hho[100];

    TCanvas *c;
    TCanvas *co;
    TLegend *leg;
    TLegend *lego = DrawMyLeg();
    TLatex *ll;

    int NEvents = tree->GetEntries();
    int IDIndex = -1;

    vector<int> TOTIDIndex;
    TOTIDIndex.clear();
    for (int j = 0; j < N; j++)
    {

        bool LoopOn = 1;
        hLowt->Reset();
        hHight->Reset();
        hRise->Reset();
        hLowU->Reset();
        hHighU->Reset();
        hLowAT->Reset();
        hHighAT->Reset();
        hLowRT->Reset();
        hHighRT->Reset();
        //h2dU->Reset();

        hLowt_Pick->Reset();
        hHight_Pick->Reset();
        hRise_Pick->Reset();
        hLowU_Pick->Reset();
        hHighU_Pick->Reset();

        tree->GetEntry(0);

        for (int i = 0; i < NEvents; i++)
        {
            //for(int i=0; i<100; i++){
            tree->GetEntry(i);
            if (j == 0)
            {
                double PMTtime = 0;
                double PMTtime_Pick = 0;
                int PMTcounter = 0;
                int PMTcounter_Pick = 0;
                for (int k = 0; k < N; k++)
                {
                    if (ChID == T0ELEchid[k] && LowTOT > UcutL)
                    {
                        PMTtime += LowThTime;
                        PMTcounter++;
                        hPMTID->Fill(T0Detchid[k]);
                        if (LowTOT < UcutR)
                        {
                            PMTcounter_Pick++;
                            PMTtime_Pick += LowThTime;
                            hPMTID_Pick->Fill(T0Detchid[k]);
                        }
                    }
                }
                hNPMT->Fill(PMTcounter);
                h2dPMT->Fill(PMTcounter, PMTtime / PMTcounter);
                hNPMT_Pick->Fill(PMTcounter_Pick);
                h2dPMT_Pick->Fill(PMTcounter_Pick, PMTtime_Pick / PMTcounter_Pick);
            }

            //else
            //LoopSwitch = 0;
            double LowThTimediff;
            double HighThTimediff;
            double risetime;
            if (ChID == T0ELEchid[j] && LowTOT > UcutL)
            {
                LowThTimediff = LowThTime;
                HighThTimediff = HighThTime;
                risetime = HighThTime - LowThTime;
                hLowt->Fill(LowThTimediff);
                hHight->Fill(HighThTimediff);
                hRise->Fill(risetime);
                hLowU->Fill(LowTOT);
                hHighU->Fill(HighTOT);
                hLowAT->Fill(LowTOT, LowThTimediff);
                hHighAT->Fill(HighTOT, LowThTimediff);
                hLowRT->Fill(LowTOT, risetime);
                hHighRT->Fill(HighTOT, risetime);
                h2dU->Fill(LowTOT, HighTOT);

                //if (LowTOT[T0ELEchid[j]] < UcutR && LowTOT[T0ELEchid[j]] > UcutL && (LowTOT[T0ELEchid[j]] > LowTOT[other[0]] || LowTOT[T0ELEchid[j]] > LowTOT[other[1]] || LowTOT[T0ELEchid[j]] > LowTOT[other[2]]))
                if (LowTOT < UcutR)

                {
                    hLowt_Pick->Fill(LowThTimediff);
                    hHight_Pick->Fill(HighThTimediff);
                    hRise_Pick->Fill(risetime);
                    hLowU_Pick->Fill(LowTOT);
                    hHighU_Pick->Fill(HighTOT);
                }
            }
        }
        cout << "check !" << endl;
        TString PngName = fileName.Copy();
        PngName = PngName.Remove(PngName.Length() - 5, 5) + "T" + T0Detchid[j];
        TString cName = Form("T%d", T0Detchid[j]);

        int counter = 0;
        sprintf(buff, "T%d", T0Detchid[j]);
        ll = DrawMyLatex(buff);

        //
        // ---------draw LowTOT vs HighTOT --------//
        //
        c = cdC(counter++);
        //ht->Rebin(2);
        h2dU->Draw("colz");
        h2dU->GetXaxis()->SetNdivisions(505);
        DrawMy2dHist(h2dU, "", "", 1, 2);
        //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
        //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
        //la = DrawMyLatex(buff, 0.2, 0.4);
        ll->Draw();
        //sprintf(buff, "%s2dU.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_2dU", cName.Data()));
        //return;

        //
        // *** Low Threshold Time distribution ***
        // **
        c = cdC(counter++);
        hLowt->Draw();
        DrawMyHist(hLowt, "", "", 1, 2, 1);
        DrawMyHist(hLowt_Pick, "", "", 2, 2, 7);
        hLowt_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hLowt, "No cut", "l");
        leg->AddEntry(hLowt_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sLowThTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowThTime", cName.Data()));

        //
        // *** HighThreshold Time distribution ***
        // **
        c = cdC(counter++);
        //c[0]->SetLogy();
        hHight->Draw();
        DrawMyHist(hHight, "", "", 1, 2, 1);
        DrawMyHist(hHight_Pick, "", "", 2, 2, 7);
        hHight_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hHight, "No cut", "l");
        leg->AddEntry(hHight_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sHighThTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighThTime", cName.Data()));

        //
        // *** Rise Time distribution ***
        // *********
        c = cdC(counter++);
        hRise->Draw();
        DrawMyHist(hRise, "", "", 1, 2, 1);
        DrawMyHist(hRise_Pick, "", "", 2, 2, 7);
        hRise_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hRise, "No cut", "l");
        leg->AddEntry(hRise_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();

        //sprintf(buff, "%sRiseTime.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_RiseTime", cName.Data()));

        //
        // *** Low TOT distribution ***
        // *********
        c = cdC(counter++);
        hLowU->GetXaxis()->SetRangeUser(UL, UR);
        hLowU->Draw();

        DrawMyHist(hLowU, "", "", 1, 2, 1);
        DrawMyHist(hLowU_Pick, "", "", 2, 2, 7);
        hLowU_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hLowU, "No cut", "l");
        leg->AddEntry(hLowU_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sLowTOTL.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowTOTL", cName.Data()));

        hLowU->GetXaxis()->SetRangeUser(0, 1.5e3);
        //sprintf(buff, "%sLowTOTS.png", PngName.Data());
        //c->SaveAs(buff);
        //opfile->WriteObject(c, Form("%s_LowTOTS", cName.Data()));

        sprintf(buff, "ho%d", T0Detchid[j]);
        TOTIDIndex.push_back(T0Detchid[j]);
        ho[T0Detchid[j]] = (TH1D *)hLowU->Clone(buff);

        //
        // *** High TOT distribution ***
        // *********
        c = cdC(counter++);
        hHighU->GetXaxis()->SetRangeUser(UL, UR);
        hHighU->Draw();

        DrawMyHist(hHighU, "", "", 1, 2, 1);
        DrawMyHist(hHighU_Pick, "", "", 2, 2, 7);
        hHighU_Pick->Draw("SAME");
        leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
        leg->AddEntry(hHighU, "No cut", "l");
        leg->AddEntry(hHighU_Pick, "With cut", "l");
        leg->Draw();

        ll->Draw();
        //sprintf(buff, "%sHighTOTL.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighTOTL", cName.Data()));

        hHighU->GetXaxis()->SetRangeUser(0, 1.5e3);
        //sprintf(buff, "%sHighTOTS.png", PngName.Data());
        //c->SaveAs(buff);

        sprintf(buff, "hho%d", T0Detchid[j]);
        hho[T0Detchid[j]] = (TH1D *)hHighU->Clone(buff);

        //
        // *** LowTOT distribution - Risetime ***
        // *********
        c = cdC(counter++);
        hLowRT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hLowRT, "", "");

        ll->Draw();
        //sprintf(buff, "%sLowRT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowRT", cName.Data()));

        //
        // *** HighTOT distribution - Risetime ***
        // *********
        c = cdC(counter++);
        hHighRT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hHighRT, "", "");

        ll->Draw();
        //sprintf(buff, "%sHighRT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighRT", cName.Data()));

        //
        // *** LowThTime-LowTOT distribution ***
        // *********
        c = cdC(counter++);
        hLowAT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
        DrawMy2dHist(hLowAT, "", "");

        ll->Draw();
        //sprintf(buff, "%sLowAT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_LowAT", cName.Data()));

        //
        // *** LowThTime-HighTOT distribution ***
        // *********
        c = cdC(counter++);
        hHighAT->Draw("colz");
        SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);

        DrawMy2dHist(hHighAT, "", "");

        ll->Draw();
        //sprintf(buff, "%sHighAT.png", PngName.Data());
        //c->SaveAs(buff);
        opfile->WriteObject(c, Form("%s_HighAT", cName.Data()));

        cout << "loop state:" << LoopOn << ", No. " << j << endl;
        if (!LoopOn)
            break;
    }

    TString PngName = fileName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5);
    //
    // ---------draw tot of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 3);
    TCanvas *p1;
    TCanvas *p2;
    TCanvas *p3;
    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowTOTL", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighTOTL", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_2dU", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
    }
    sprintf(buff, "%sTOT_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw time of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 3);

    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowThTime", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighThTime", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_RiseTime", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
    }
    sprintf(buff, "%sTime_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw 2d relationship  of PMT--------//
    //

    c = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1000);
    c->Divide(N, 4);

    TCanvas *p4;
    for (int i = 0; i < N; i++)
    {
        p1 = (TCanvas *)opfile->Get(Form("T%d_LowRT", T0Detchid[i]));
        p2 = (TCanvas *)opfile->Get(Form("T%d_HighRT", T0Detchid[i]));
        p3 = (TCanvas *)opfile->Get(Form("T%d_LowAT", T0Detchid[i]));
        p4 = (TCanvas *)opfile->Get(Form("T%d_HighAT", T0Detchid[i]));

        c->cd(T0Detchid[i] + 1);
        p1->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 1 * N);
        p2->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 2 * N);
        p3->DrawClonePad();
        c->cd(T0Detchid[i] + 1 + 3 * N);
        p4->DrawClonePad();
    }
    sprintf(buff, "%s2Drelationship_of_T0.png", PngName.Data());
    c->SaveAs(buff);

    //
    // ---------draw efficiency of PMT--------//
    //
    double PMTEff[N];
    double PMTEff_Pick[N];
    c = cdC(counter++);
    DrawMyHist(hPMTID, "", "", 1, 2);
    //ht->Rebin(2);
    hPMTID->Draw();
    hPMTID->GetXaxis()->SetNdivisions(505);
    DrawMyHist(hPMTID_Pick, "", "", 2, 2, 7);
    hPMTID_Pick->Draw("same");
    leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(hPMTID, "No cut", "l");
    leg->AddEntry(hPMTID_Pick, "With cut", "l");
    leg->Draw();
    for (int i = 0; i < N; i++)
    {

        PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / NEvents;
        PMTEff_Pick[i] = hPMTID_Pick->GetBinContent(hPMTID_Pick->FindBin(i)) / NEvents;
        sprintf(buff, "PMTID=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i, PMTEff[i] * 100, PMTEff_Pick[i] * 100);
        ll = DrawMyLatex(buff, 0.3, 0.4 + 0.1 * i, 42, 0.06);
        ll->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_PMT.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_PMT"));

    //return;

    //
    // *** Number of fired PMT ***
    // **
    double Eff[N + 1];
    double Eff_Pick[N + 1];
    TLatex *la;

    c = cdC(counter++);
    DrawMyHist(hNPMT, "", "", 1, 2);
    //ht->Rebin(2);
    hNPMT->Draw();
    hNPMT->GetXaxis()->SetNdivisions(505);
    DrawMyHist(hNPMT_Pick, "", "", 2, 2, 7);
    hNPMT_Pick->Draw("SAME");
    leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(hNPMT, "No cut", "l");
    leg->AddEntry(hNPMT_Pick, "With cut", "l");
    leg->Draw();
    for (int i = 0; i < N + 1; i++)
    {
        Eff[i] = hNPMT->GetBinContent(hNPMT->FindBin(i)) / hNPMT->Integral();
        Eff_Pick[i] = hNPMT_Pick->GetBinContent(hNPMT_Pick->FindBin(i)) / hNPMT_Pick->Integral();
        sprintf(buff, "NPMT=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i, Eff[i] * 100, Eff_Pick[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.4 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_T0.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_T0"));
    //
    // ---------draw NPMT vs TR --------//
    //
    c = cdC(counter++);
    DrawMy2dHist(h2dPMT, "", "", 1, 2);
    //ht->Rebin(2);
    h2dPMT->Draw("colz");
    h2dPMT->GetXaxis()->SetNdivisions(505);
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    //sprintf(buff, "%s2dPMT.png", PngName.Data());
    //c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_T0_vs_Time"));
    //return;

    int temp;
    int max_y_content;
    sort(TOTIDIndex.begin(), TOTIDIndex.end());

    //
    // ---------draw LowTOT together --------//
    //
    co = cdC(counter++);
    temp = 0;
    max_y_content = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        FillHistArea(ho[TOTIDIndex.at(i)], clr[i], 0.8, fst[2]);
        sprintf(buff, "T%d", TOTIDIndex.at(i));
        lego->AddEntry(ho[TOTIDIndex.at(i)], buff, "lpf");
        temp = ho[TOTIDIndex.at(i)]->GetMaximum();
        max_y_content = max_y_content > temp ? max_y_content : temp;
        if (i == 0)
            ho[TOTIDIndex.at(i)]->Draw();
        else
            ho[TOTIDIndex.at(i)]->Draw("same");
    }
    ho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max_y_content * 1.1);
    lego->Draw();
    sprintf(buff, "%sLowTOTDrawtogether.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_lowtot"));

    //
    // ---------draw HighTOT together --------//
    //
    lego->Clear();
    co = cdC(counter++);
    temp = 0;
    max_y_content = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        FillHistArea(hho[TOTIDIndex.at(i)], clr[i], 0.1, fst[2]);
        sprintf(buff, "T%d", TOTIDIndex.at(i));
        lego->AddEntry(hho[TOTIDIndex.at(i)], buff, "lpf");
        temp = hho[TOTIDIndex.at(i)]->GetMaximum();
        max_y_content = max_y_content > temp ? max_y_content : temp;

        if (i == 0)
            hho[TOTIDIndex.at(i)]->Draw();
        else
            hho[TOTIDIndex.at(i)]->Draw("same");
    }
    hho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max_y_content * 1.1);
    lego->Draw();
    sprintf(buff, "%sHighTOTDrawtogether.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_hightot"));
}

void PrintRawFTOF(TString fileName, TString opfileName, int force = 0)
{

    //cout<<"FTOFChNum: "<<FTOFChNum<<endl;
    //return;

    TFile *opfile;
    FileStat_t fStat;
    gSystem->GetPathInfo(opfileName, fStat);
    if (fStat.fSize != 0 && force != 1)
    {

        cout << "ERROR! " << opfileName << " is exist!" << endl;
        opfile = new TFile(opfileName, "read");
    }
    else
        opfile = new TFile(opfileName, "recreate");
    cout << "====>>  Create the root file : " << opfileName << endl;

    TGaxis::SetMaxDigits(3);

    Color_t clr[] = {1, kRed + 3, kRed, kOrange + 2, kYellow - 3, kSpring - 1, kGreen + 3, kCyan + 1, kBlue, kViolet + 2, kMagenta - 4, kPink - 9, kGray + 2}; //color
    Style_t fst[] = {1001, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3144, 3244, 3481, 3315, 3351};                                                            //

    //sprintf(path, "/mnt/f/CRsys/CosmicRayTest/%s%s", date, root);
    //sprintf(path, "/mnt/d/ExpDATA/labtest/CosmicRay/%s%s", date, root);

    //sprintf(rootname,"%s/beamtest%s/Analysismatch%s",rootpath,date,date);

    //double LowThTimediff;
    double tL = -800e3, tR = -600e3;
    double rtL = -1e3, rtR = 0.2e3;
    double TcutL = 1;
    double TcutR = 2e3; //reffered Tube;
    //TcutL=1;
    //TcutR=4e3; //reffered Tube;

    double UL = 0, UR = 30e3;
    double UcutL = 1, UcutR = 3e3; //reffered Tube;

    //UcutL=1;
    //UcutR=4e3; //reffered Tube;

    //int rbt=10;
    //int rbU=8;
    TH1D *h1d;
    TH1D *h1dpick;
    TH1D *h2d;
    TH1D *h2dpick;
    TFile *file = new TFile(fileName, "read");
    //TTree* tree =  (TTree*)file->Get("Analysis");
    TTree *tree = (TTree *)file->Get("tree");
    cout << "====>>  Start to open the root file : " << fileName << endl;

    int TrackID;
    int ChID;
    double LowThTime;
    double HighThTime;
    double LowTOT;
    double HighTOT;
    tree->SetBranchAddress("TrackID", &TrackID);
    tree->SetBranchAddress("ChID", &ChID);
    tree->SetBranchAddress("LowThTime", &LowThTime);
    tree->SetBranchAddress("HighThTime", &HighThTime);
    tree->SetBranchAddress("LowTOT", &LowTOT);
    tree->SetBranchAddress("HighTOT", &HighTOT);

    TH1D *hLowt = new TH1D("hLowt", "Time Resolution;LowTh time (ps);Counts", 500, tL, tR);
    TH1D *hLowt_Pick = (TH1D *)hLowt->Clone("hLowt_Pick");

    TH1D *hHight = new TH1D("hHight", "Time Resolution;HighTh time (ps);Counts", 500, tL, tR);
    TH1D *hHight_Pick = (TH1D *)hHight->Clone("hHight_Pick");

    TH1D *hRise = new TH1D("hRise", "Rise Time; Rise time (ps);Counts", 0.2e3, rtL, rtR);
    TH1D *hRise_Pick = (TH1D *)hRise->Clone("hRise_Pick");

    TH1D *hHighU = new TH1D("hHighU", "HighTOT;HighTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hLowU = new TH1D("hLowU", "LowTOT;LowTOT (ps);Counts", 5e3, UL, UR);
    TH1D *hHighU_Pick = (TH1D *)hHighU->Clone("hHighU_Pick");
    TH1D *hLowU_Pick = (TH1D *)hLowU->Clone("hLowU_Pick");
    TH2D *h2dU = new TH2D("h2dU", ";LowTOT (ps)T;HighTOT (ps)", 500, UL, UR, 500, UL, UR);

    TH1D *hPMTID = new TH1D("hPMTID", ";PMTID;Counts", FTOFChNum, 0, FTOFChNum);
    TH1D *hPMTID_Pick = (TH1D *)hPMTID->Clone("hPMTID_Pick");

    TH1D *hNPMT = new TH1D("hNPMT", ";Number of Fired PMT;Counts", FTOFChNum + 1, 0, FTOFChNum + 1);
    TH1D *hNPMT_Pick = (TH1D *)hNPMT->Clone("hNPMT_Pick");

    TH2D *h2dPMT = new TH2D("h2dPMT", ";Number of Fired PMT;Timediff (ps)", FTOFChNum + 1, 0, FTOFChNum + 1, 100, tL, tR);
    TH2D *h2dPMT_Pick = (TH2D *)h2dPMT->Clone("h2dPMT_Pick");
    TH2D *hHighAT = new TH2D(buff, "hHighAT;HighTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowAT = new TH2D(buff, "hLowAT;LowTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowRT = new TH2D(buff, "hLowRT;LowTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);
    TH2D *hHighRT = new TH2D(buff, "hHighRT;HighTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);

    TH1D *ho[100];
    TH1D *hho[100];

    TCanvas *c;
    TCanvas *co;
    TLegend *leg;
    TLegend *lego = DrawMyLeg();
    TLatex *ll;

    int NEvents = tree->GetEntries();
    int IDIndex = -1;
    vector<int> TOTIDIndex;
    TOTIDIndex.clear();
    for (int j = 0; j < FTOFChNum; j++)
    {
        if (Eabled[j])
        {

            bool LoopOn = 1;
            hLowt->Reset();
            hHight->Reset();
            hRise->Reset();
            hLowU->Reset();
            hHighU->Reset();
            hLowAT->Reset();
            hHighAT->Reset();
            hLowRT->Reset();
            hHighRT->Reset();
            //h2dU->Reset();

            hLowt_Pick->Reset();
            hHight_Pick->Reset();
            hRise_Pick->Reset();
            hLowU_Pick->Reset();
            hHighU_Pick->Reset();

            tree->GetEntry(0);

            for (int i = 0; i < NEvents; i++)
            {
                //for(int i=0; i<100; i++){
                tree->GetEntry(i);
                if (j == 0)
                {
                    double PMTtime = 0;
                    double PMTtime_Pick = 0;
                    int PMTcounter = 0;
                    int PMTcounter_Pick = 0;
                    for (int k = 0; k < FTOFChNum; k++)
                    {
                        if (ChID == FTOFELEchid[k] && Eabled[k])
                        {

                            if (LowTOT > UcutL)
                            {
                                PMTtime += LowThTime;
                                PMTcounter++;
                                hPMTID->Fill(FTOFDetchid[k]);
                                if (LowTOT < UcutR)
                                {
                                    PMTcounter_Pick++;
                                    PMTtime_Pick += LowThTime;
                                    hPMTID_Pick->Fill(FTOFDetchid[k]);
                                }
                            }
                        }
                    }
                    hNPMT->Fill(PMTcounter);
                    h2dPMT->Fill(PMTcounter, PMTtime / PMTcounter);
                    hNPMT_Pick->Fill(PMTcounter_Pick);
                    h2dPMT_Pick->Fill(PMTcounter_Pick, PMTtime_Pick / PMTcounter_Pick);
                }

                //else
                //LoopSwitch = 0;
                double LowThTimediff;
                double HighThTimediff;
                double risetime;
                if (ChID == FTOFELEchid[j] && LowTOT > UcutL)
                {
                    LowThTimediff = LowThTime;
                    HighThTimediff = HighThTime;
                    risetime = HighThTime - LowThTime;
                    hLowt->Fill(LowThTimediff);
                    hHight->Fill(HighThTimediff);
                    hRise->Fill(risetime);
                    hLowU->Fill(LowTOT);
                    hHighU->Fill(HighTOT);
                    hLowAT->Fill(LowTOT, LowThTimediff);
                    hHighAT->Fill(HighTOT, LowThTimediff);
                    hLowRT->Fill(LowTOT, risetime);
                    hHighRT->Fill(HighTOT, risetime);
                    h2dU->Fill(LowTOT, HighTOT);

                    //if (LowTOT[j] < UcutR && LowTOT[j] > UcutL && (LowTOT[j] > LowTOT[other[0]] || LowTOT[j] > LowTOT[other[1]] || LowTOT[j] > LowTOT[other[2]]))
                    if (LowTOT < UcutR)

                    {
                        hLowt_Pick->Fill(LowThTimediff);
                        hHight_Pick->Fill(HighThTimediff);
                        hRise_Pick->Fill(risetime);
                        hLowU_Pick->Fill(LowTOT);
                        hHighU_Pick->Fill(HighTOT);
                    }
                }
            }
            cout << "check !" << endl;
            TString PngName = fileName.Copy();
            PngName = PngName.Remove(PngName.Length() - 5, 5) + "CH" + FTOFDetchid[j];
            TString cName = Form("CH%d", FTOFDetchid[j]);

            int counter = 0;
            sprintf(buff, "Ch%d", FTOFDetchid[j]);
            ll = DrawMyLatex(buff);

            //
            // ---------draw LowTOT vs HighTOT --------//
            //
            c = cdC(counter++);
            //ht->Rebin(2);
            h2dU->Draw("colz");
            h2dU->GetXaxis()->SetNdivisions(505);
            DrawMy2dHist(h2dU, "", "", 1, 2);
            //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
            //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
            //la = DrawMyLatex(buff, 0.2, 0.4);
            ll->Draw();
            //sprintf(buff, "%s2dU.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_2dU", cName.Data()));
            //return;

            //
            // *** Low Threshold Time distribution ***
            // **
            c = cdC(counter++);
            hLowt->Draw();
            DrawMyHist(hLowt, "", "", 1, 2, 1);
            DrawMyHist(hLowt_Pick, "", "", 2, 2, 7);
            hLowt_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hLowt, "No cut", "l");
            leg->AddEntry(hLowt_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sLowThTime.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_LowThTime", cName.Data()));

            //
            // *** HighThreshold Time distribution ***
            // **
            c = cdC(counter++);
            //c[0]->SetLogy();
            hHight->Draw();
            DrawMyHist(hHight, "", "", 1, 2, 1);
            DrawMyHist(hHight_Pick, "", "", 2, 2, 7);
            hHight_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hHight, "No cut", "l");
            leg->AddEntry(hHight_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sHighThTime.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_HighThTime", cName.Data()));

            //
            // *** Rise Time distribution ***
            // *********
            c = cdC(counter++);
            hRise->Draw();
            DrawMyHist(hRise, "", "", 1, 2, 1);
            DrawMyHist(hRise_Pick, "", "", 2, 2, 7);
            hRise_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hRise, "No cut", "l");
            leg->AddEntry(hRise_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();

            //sprintf(buff, "%sRiseTime.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_RiseTime", cName.Data()));

            //
            // *** Low TOT distribution ***
            // *********
            c = cdC(counter++);
            hLowU->GetXaxis()->SetRangeUser(UL, UR);
            hLowU->Draw();

            DrawMyHist(hLowU, "", "", 1, 2, 1);
            DrawMyHist(hLowU_Pick, "", "", 2, 2, 7);
            hLowU_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hLowU, "No cut", "l");
            leg->AddEntry(hLowU_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sLowTOTL.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_LowTOTL", cName.Data()));

            hLowU->GetXaxis()->SetRangeUser(0, 1.5e3);
            //sprintf(buff, "%sLowTOTS.png", PngName.Data());
            //c->SaveAs(buff);
            //opfile->WriteObject(c, Form("%s_LowTOTS", cName.Data()));

            sprintf(buff, "ho%d", FTOFDetchid[j]);
            TOTIDIndex.push_back(FTOFDetchid[j]);
            ho[FTOFDetchid[j]] = (TH1D *)hLowU->Clone(buff);
            //
            // *** High TOT distribution ***
            // *********
            c = cdC(counter++);
            hHighU->GetXaxis()->SetRangeUser(UL, UR);
            hHighU->Draw();

            DrawMyHist(hHighU, "", "", 1, 2, 1);
            DrawMyHist(hHighU_Pick, "", "", 2, 2, 7);
            hHighU_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hHighU, "No cut", "l");
            leg->AddEntry(hHighU_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sHighTOTL.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_HighTOTL", cName.Data()));

            hHighU->GetXaxis()->SetRangeUser(0, 1.5e3);
            //sprintf(buff, "%sHighTOTS.png", PngName.Data());
            //c->SaveAs(buff);

            sprintf(buff, "hho%d", FTOFDetchid[j]);
            hho[FTOFDetchid[j]] = (TH1D *)hHighU->Clone(buff);

            //
            // *** LowTOT distribution - Risetime ***
            // *********
            c = cdC(counter++);
            hLowRT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hLowRT, "", "");

            ll->Draw();
            //sprintf(buff, "%sLowRT.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_LowRT", cName.Data()));

            //
            // *** HighTOT distribution - Risetime ***
            // *********
            c = cdC(counter++);
            hHighRT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hHighRT, "", "");

            ll->Draw();
            //sprintf(buff, "%sHighRT.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_HighRT", cName.Data()));

            //
            // *** LowThTime-LowTOT distribution ***
            // *********
            c = cdC(counter++);
            hLowAT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hLowAT, "", "");

            ll->Draw();
            //sprintf(buff, "%sLowAT.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_LowAT", cName.Data()));

            //
            // *** LowThTime-HighTOT distribution ***
            // *********
            c = cdC(counter++);
            hHighAT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);

            DrawMy2dHist(hHighAT, "", "");

            ll->Draw();
            //sprintf(buff, "%sHighAT.png", PngName.Data());
            //c->SaveAs(buff);
            opfile->WriteObject(c, Form("%s_HighAT", cName.Data()));

            cout << "loop state:" << LoopOn << ", No. " << j << endl;
            if (!LoopOn)
                break;
        }
    }

    TString PngName = fileName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5);

    //
    // ---------draw tot of PMT--------//
    //
    TCanvas *c_all1;
    TCanvas *c_all2;
    TCanvas *c_all3;
    TCanvas *c_all4;
    TCanvas *p1;
    TCanvas *p2;
    TCanvas *p3;
    TCanvas *p4;

    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);

    c_all1->Divide(4, 4);
    c_all2->Divide(4, 4);
    c_all3->Divide(4, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowTOTL", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighTOTL", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_2dU", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowTOT_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighTOT_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%s2dTOT_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);

    //
    // ---------draw time of PMT--------//
    //
    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all1->Divide(4, 4);
    c_all2->Divide(4, 4);
    c_all3->Divide(4, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowThTime", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighThTime", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_RiseTime", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowThtime_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighThtime_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%sRiseTime_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);

    //
    // ---------draw 2d relationship  of PMT--------//
    //
    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all4 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all1->Divide(4, 4);
    c_all2->Divide(4, 4);
    c_all3->Divide(4, 4);
    c_all4->Divide(4, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowRT", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighRT", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_LowAT", FTOFDetchid[i]));
            p4 = (TCanvas *)opfile->Get(Form("CH%d_HighAT", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
            c_all4->cd(FTOFDetchid[i] + 1);
            p4->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowRT_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighRT_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%sLowAT_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);
    sprintf(buff, "%sHighAT_of_FTOF.png", PngName.Data());
    c_all4->SaveAs(buff);

    //
    // ---------draw efficiency of PMT--------//
    //
    double PMTEff[FTOFChNum];
    double PMTEff_Pick[FTOFChNum];
    c = cdC(counter++);
    DrawMyHist(hPMTID, "", "", 1, 2);
    //ht->Rebin(2);
    hPMTID->Draw();
    hPMTID->GetXaxis()->SetNdivisions(505);
    DrawMyHist(hPMTID_Pick, "", "", 2, 2, 7);
    hPMTID_Pick->Draw("same");
    leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(hPMTID, "No cut", "l");
    leg->AddEntry(hPMTID_Pick, "With cut", "l");
    leg->Draw();
    for (int i = 0; i < FTOFChNum; i++)
    {

        PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / NEvents;
        PMTEff_Pick[i] = hPMTID_Pick->GetBinContent(hPMTID_Pick->FindBin(i)) / NEvents;
        sprintf(buff, "PMTID=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i + 1, PMTEff[i] * 100, PMTEff_Pick[i] * 100);
        ll = DrawMyLatex(buff, 0.3, 0.4 + 0.1 * i, 42, 0.06);
        ll->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_FTOFCH.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_PMT"));

    //return;

    //
    // *** Number of fired PMT ***
    // **
    double Eff[FTOFChNum + 1];
    double Eff_Pick[FTOFChNum + 1];
    TLatex *la;

    c = cdC(counter++);
    DrawMyHist(hNPMT, "", "", 1, 2);
    //ht->Rebin(2);
    hNPMT->Draw();
    hNPMT->GetXaxis()->SetNdivisions(505);
    DrawMyHist(hNPMT_Pick, "", "", 2, 2, 7);
    hNPMT_Pick->Draw("SAME");
    leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    leg->AddEntry(hNPMT, "No cut", "l");
    leg->AddEntry(hNPMT_Pick, "With cut", "l");
    leg->Draw();
    for (int i = 0; i < FTOFChNum + 1; i++)
    {
        Eff[i] = hNPMT->GetBinContent(hNPMT->FindBin(i)) / hNPMT->Integral();
        Eff_Pick[i] = hNPMT_Pick->GetBinContent(hNPMT_Pick->FindBin(i)) / hNPMT_Pick->Integral();
        sprintf(buff, "NPMT=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i, Eff[i] * 100, Eff_Pick[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.4 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_FTOF.png", PngName.Data());
    c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_FTOF"));
    //
    // ---------draw NPMT vs TR --------//
    //
    c = cdC(counter++);
    DrawMy2dHist(h2dPMT, "", "", 1, 2);
    //ht->Rebin(2);
    h2dPMT->Draw("colz");
    h2dPMT->GetXaxis()->SetNdivisions(505);
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    //sprintf(buff, "%s2dPMT.png", PngName.Data());
    //c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_FTOF_vs_Time"));
    //return;

    //
    // ---------draw LowTOT together --------//
    //
    co = cdC(counter++);
    int max = 0;
    int temp = 0;
    sort(TOTIDIndex.begin(), TOTIDIndex.end());
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        cout << "lowTOThist: " << i << endl;
        temp = ho[TOTIDIndex.at(i)]->GetMaximum();
        max = max < temp ? temp : max;
        FillHistArea(ho[TOTIDIndex.at(i)], clr[i], 0.8, fst[2]);
        sprintf(buff, "Ch%d", TOTIDIndex.at(i));
        lego->AddEntry(ho[TOTIDIndex.at(i)], buff, "lpf");
        if (i == 0)
            ho[TOTIDIndex.at(i)]->Draw();

        else
            ho[TOTIDIndex.at(i)]->Draw("same");
    }
    ho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max * 1.1);
    lego->Draw();
    sprintf(buff, "%sLowTOTDrawtogether_of_FTOF.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_lowtot"));

    //
    // ---------draw HighTOT together --------//
    //
    lego->Clear();
    co = cdC(counter++);
    max = 0;
    temp = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {

        cout << "highTOThist: " << i << endl;
        temp = hho[TOTIDIndex.at(i)]->GetMaximum();
        max = max < temp ? temp : max;
        FillHistArea(hho[TOTIDIndex.at(i)], clr[i], 0.1, fst[2]);
        sprintf(buff, "Ch%d", TOTIDIndex.at(i));
        lego->AddEntry(hho[TOTIDIndex.at(i)], buff, "lpf");
        if (i == 0)
            hho[TOTIDIndex.at(i)]->Draw();
        else
            hho[TOTIDIndex.at(i)]->Draw("same");
    }
    hho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max * 1.1);
    lego->Draw();
    sprintf(buff, "%sHighTOTDrawtogether__of_FTOF.png", PngName.Data());
    co->SaveAs(buff);
    opfile->WriteObject(c, Form("all_ch_hightot"));
}
*/
void PrintFTOF(TString fileName, TString opfileName, int force = 0)
{
    cout << "FTOFChNum: " << FTOFChNum << endl;
    //return;
    //ReadFTOFmap("FTOFdet2ELE_2ndmap.csv");
    //ReadFTOFmap("FTOFdet2ELE.csv");
    TFile *opfile;
    FileStat_t fStat;
    bool existflag = 0;
    gSystem->GetPathInfo(opfileName, fStat);
    if (fStat.fSize != 0 && force != 1)
    {

        cout << "ERROR! " << opfileName << " is exist!" << endl;
        opfile = new TFile(opfileName, "read");
        existflag = 1;
    }
    else
        opfile = new TFile(opfileName, "recreate");
    cout << "====>>  Create the root file : " << opfileName << endl;

    TGaxis::SetMaxDigits(3);

    Color_t clr[] = {1, kRed + 3, kRed, kOrange + 2, kYellow - 3, kSpring - 1, kGreen + 3, kCyan + 1, kBlue, kViolet + 2, kMagenta - 4, kPink - 9, kGray + 2}; //color
    Style_t fst[] = {1001, 3001, 3002, 3003, 3004, 3005, 3006, 3007, 3144, 3244, 3481, 3315, 3351};                                                            //

    //sprintf(path, "/mnt/f/CRsys/CosmicRayTest/%s%s", date, root);
    //sprintf(path, "/mnt/d/ExpDATA/labtest/CosmicRay/%s%s", date, root);

    //sprintf(rootname,"%s/beamtest%s/Analysismatch%s",rootpath,date,date);

    //double LowThTimediff;
    double tL = -500e3, tR = -400e3;
    double rtL = -1e3, rtR = 0.2e3;
    double TcutL = 1;
    double TcutR = 2e3; //reffered Tube;
    //TcutL=1;
    //TcutR=4e3; //reffered Tube;

    double UL = 0, UR = 5e3;
    double UcutL = 1, UcutR = 3e3; //reffered Tube;

    //UcutL=1;
    //UcutR=4e3; //reffered Tube;
    for (int i = 0; i < FTOFChNum; i++)
    {
        FTOFELEchid[i] = i;
        FTOFDetchid[i] = i;
        Eabled[i] = 1;
    }
    //int rbt=10;
    //int rbU=8;
    TH1D *h1d;
    TH1D *h1dpick;
    TH1D *h2d;
    TH1D *h2dpick;
    TFile *file = new TFile(fileName, "read");
    //TTree* tree =  (TTree*)file->Get("Analysis");
    TTree *tree = (TTree *)file->Get("tree");
    cout << "====>>  Start to open the root file : " << fileName << endl;

    FTOFHit *FTOFhit = new FTOFHit;
    tree->SetBranchAddress("FTOFhit", &FTOFhit);

    TH1D *hLowt = new TH1D("hLowt", "Time Resolution;LowTh time (ps);Counts", 500, tL, tR);
    TH1D *hLowt_Pick = (TH1D *)hLowt->Clone("hLowt_Pick");

    TH1D *hHight = new TH1D("hHight", "Time Resolution;HighTh time (ps);Counts", 500, tL, tR);
    TH1D *hHight_Pick = (TH1D *)hHight->Clone("hHight_Pick");

    TH1D *hRise = new TH1D("hRise", "Rise Time; Rise time (ps);Counts", 0.2e3, rtL, rtR);
    TH1D *hRise_Pick = (TH1D *)hRise->Clone("hRise_Pick");

    TH1D *hHighU = new TH1D("hHighU", "HighTOT;HighTOT (ps);Counts", 1e3, UL, UR);
    TH1D *hLowU = new TH1D("hLowU", "LowTOT;LowTOT (ps);Counts", 1e3, UL, UR);
    TH1D *hHighU_Pick = (TH1D *)hHighU->Clone("hHighU_Pick");
    TH1D *hLowU_Pick = (TH1D *)hLowU->Clone("hLowU_Pick");
    TH2D *h2dU = new TH2D("h2dU", ";LowTOT (ps)T;HighTOT (ps)", 500, UL, UR, 500, UL, UR);

    TH1D *hPMTID = new TH1D("hPMTID", ";PMTID;Counts", FTOFChNum, 0, FTOFChNum);
    TH1D *hPMTID_Pick = (TH1D *)hPMTID->Clone("hPMTID_Pick");

    TH1D *hNPMT = new TH1D("hNPMT", ";Number of Fired PMT;Counts", FTOFChNum + 1, 0, FTOFChNum + 1);
    TH1D *hNPMT_Pick = (TH1D *)hNPMT->Clone("hNPMT_Pick");

    TH2D *h2dPMT = new TH2D("h2dPMT", ";Number of Fired PMT;Timediff (ps)", FTOFChNum + 1, 0, FTOFChNum + 1, 100, tL, tR);
    TH2D *h2dPMT_Pick = (TH2D *)h2dPMT->Clone("h2dPMT_Pick");
    TH2D *hHighAT = new TH2D("", "hHighAT;HighTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowAT = new TH2D("", "hLowAT;LowTOT (ps); Timediff (ps)", 100, UcutL, 1.5e3, 100, tL, tR);
    TH2D *hLowRT = new TH2D("", "hLowRT;LowTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);
    TH2D *hHighRT = new TH2D("", "hHighRT;HighTOT (ps); Risetime (ps)", 100, UcutL, 1.5e3, 100, rtL, rtR);
    //TH2I *hpos = new TH2I("hpos", "hpos; x ;  y", 8, 0, 8, 4, 0, 4);
    TH2I *hpos = new TH2I("hpos", "hpos; x ;  y", 24, 0, 24, 4, 0, 4);

    TH1D *ho[100];
    TH1D *hho[100];

    TCanvas *c;
    TCanvas *co;
    TLegend *leg;
    TLegend *lego = DrawMyLeg();
    TLatex *ll;
    int eventcounter = 0;
    int NEvents = tree->GetEntries();
    int IDIndex = -1;
    int counter = 0;

    vector<int> TOTIDIndex;
    TOTIDIndex.clear();
    bool onetimeflag = 1;

    for (int j = 0; j < FTOFChNum; j++)
    {
        //if (Eabled[j])
        if (1)
        {

            bool LoopOn = 1;
            hLowt->Reset();
            hHight->Reset();
            hRise->Reset();
            hLowU->Reset();
            hHighU->Reset();
            hLowAT->Reset();
            hHighAT->Reset();
            hLowRT->Reset();
            hHighRT->Reset();
            h2dU->Reset();

            hLowt_Pick->Reset();
            hHight_Pick->Reset();
            hRise_Pick->Reset();
            hLowU_Pick->Reset();
            hHighU_Pick->Reset();

            tree->GetEntry(0);

            for (int i = 0; i < NEvents; i++)
            {
                //for(int i=0; i<100; i++){
                tree->GetEntry(i);
                if (onetimeflag)
                {
                    double PMTtime = 0;
                    double PMTtime_Pick = 0;
                    int PMTcounter = 0;
                    int PMTcounter_Pick = 0;
                    for (int k = 0; k < FTOFChNum; k++)
                    {
                        //if (Eabled[k])
                        if (1)
                        {

                            if (FTOFhit->lowtot[FTOFELEchid[k]] > UcutL)
                            {
                                PMTtime += FTOFhit->triggertimediff[FTOFELEchid[k]];
                                PMTcounter++;
                                hPMTID->Fill(FTOFDetchid[k]);
                                int x = -999;
                                int y = -999;
                                IndexFTOFpos(FTOFDetchid[k], x, y);
                                hpos->Fill(x, y);
                                if (FTOFhit->lowtot[FTOFELEchid[k]] < UcutR)
                                {
                                    PMTcounter_Pick++;
                                    PMTtime_Pick += FTOFhit->triggertimediff[FTOFELEchid[k]];
                                    hPMTID_Pick->Fill(FTOFDetchid[k]);
                                }
                            }
                        }
                    }
                    if (PMTcounter > 0)
                    {
                        eventcounter++;
                        h2dPMT->Fill(PMTcounter, PMTtime / PMTcounter);
                        h2dPMT_Pick->Fill(PMTcounter_Pick, PMTtime_Pick / PMTcounter_Pick);
                    }
                    hNPMT_Pick->Fill(PMTcounter_Pick);
                    hNPMT->Fill(PMTcounter);
                }

                //else
                //LoopSwitch = 0;
                double LowThTimediff;
                double HighThTimediff;
                double risetime;
                if (FTOFhit->lowtot[FTOFELEchid[j]] > UcutL)
                {
                    LowThTimediff = FTOFhit->triggertimediff[FTOFELEchid[j]];
                    HighThTimediff = FTOFhit->highthtime[FTOFELEchid[j]] + FTOFhit->triggertimediff[FTOFELEchid[j]] - FTOFhit->lowthtime[FTOFELEchid[j]];
                    risetime = FTOFhit->highthtime[FTOFELEchid[j]] - FTOFhit->lowthtime[FTOFELEchid[j]];
                    hLowt->Fill(LowThTimediff);
                    hHight->Fill(HighThTimediff);
                    hRise->Fill(risetime);
                    hLowU->Fill(FTOFhit->lowtot[FTOFELEchid[j]]);
                    hHighU->Fill(FTOFhit->hightot[FTOFELEchid[j]]);
                    hLowAT->Fill(FTOFhit->lowtot[FTOFELEchid[j]], LowThTimediff);
                    hHighAT->Fill(FTOFhit->hightot[FTOFELEchid[j]], LowThTimediff);
                    hLowRT->Fill(FTOFhit->lowtot[FTOFELEchid[j]], risetime);
                    hHighRT->Fill(FTOFhit->hightot[FTOFELEchid[j]], risetime);
                    h2dU->Fill(FTOFhit->lowtot[FTOFELEchid[j]], FTOFhit->hightot[FTOFELEchid[j]]);

                    //if (LowTOT[j] < UcutR && LowTOT[j] > UcutL && (LowTOT[j] > LowTOT[other[0]] || LowTOT[j] > LowTOT[other[1]] || LowTOT[j] > LowTOT[other[2]]))
                    if (FTOFhit->lowtot[FTOFELEchid[j]] < UcutR)

                    {
                        hLowt_Pick->Fill(LowThTimediff);
                        hHight_Pick->Fill(HighThTimediff);
                        hRise_Pick->Fill(risetime);
                        hLowU_Pick->Fill(FTOFhit->lowtot[FTOFELEchid[j]]);
                        hHighU_Pick->Fill(FTOFhit->hightot[FTOFELEchid[j]]);
                    }
                }
            }
            onetimeflag = 0;

            cout << "check !" << endl;
            TString PngName = fileName.Copy();
            PngName = PngName.Remove(PngName.Length() - 5, 5) + "CH" + TString::Itoa(FTOFDetchid[j], 10);
            TString cName = Form("CH%d", FTOFDetchid[j]);

            sprintf(buff, "Ch%d", FTOFDetchid[j]);
            ll = DrawMyLatex(buff);

            //
            // ---------draw LowTOT vs HighTOT --------//
            //
            c = cdC(1);
            //ht->Rebin(2);
            h2dU->Draw("colz");
            h2dU->GetXaxis()->SetNdivisions(505);
            DrawMy2dHist(h2dU, "", "", 1, 2);
            //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
            //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
            //la = DrawMyLatex(buff, 0.2, 0.4);
            ll->Draw();
            //sprintf(buff, "%s2dU.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_2dU", cName.Data()));
            c->Close();
            //return;

            //
            // *** Low Threshold Time distribution ***
            // **
            c = cdC(2);
            hLowt->Draw();
            DrawMyHist(hLowt, "", "", 1, 2, 1);
            DrawMyHist(hLowt_Pick, "", "", 2, 2, 7);
            hLowt_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hLowt, "No cut", "l");
            leg->AddEntry(hLowt_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sLowThTime.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_LowThTime", cName.Data()));
            c->Close();

            //
            // *** HighThreshold Time distribution ***
            // **
            c = cdC(3);
            //c[0]->SetLogy();
            hHight->Draw();
            DrawMyHist(hHight, "", "", 1, 2, 1);
            DrawMyHist(hHight_Pick, "", "", 2, 2, 7);
            hHight_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hHight, "No cut", "l");
            leg->AddEntry(hHight_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sHighThTime.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_HighThTime", cName.Data()));
            c->Close();

            //
            // *** Rise Time distribution ***
            // *********
            c = cdC(4);
            hRise->Draw();
            DrawMyHist(hRise, "", "", 1, 2, 1);
            DrawMyHist(hRise_Pick, "", "", 2, 2, 7);
            hRise_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hRise, "No cut", "l");
            leg->AddEntry(hRise_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();

            //sprintf(buff, "%sRiseTime.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_RiseTime", cName.Data()));
            c->Close();

            //
            // *** Low TOT distribution ***
            // *********
            c = cdC(5);
            hLowU->GetXaxis()->SetRangeUser(UL, UR);
            hLowU->Draw();

            DrawMyHist(hLowU, "", "", 1, 2, 1);
            DrawMyHist(hLowU_Pick, "", "", 2, 2, 7);
            hLowU_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hLowU, "No cut", "l");
            leg->AddEntry(hLowU_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sLowTOTL.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_LowTOTL", cName.Data()));
            c->Close();

            //hLowU->GetXaxis()->SetRangeUser(0, 1.5e3);
            //sprintf(buff, "%sLowTOTS.png", PngName.Data());
            //c->SaveAs(buff);
            //opfile->WriteObject(c, Form("%s_LowTOTS", cName.Data()));

            sprintf(buff, "ho%d", FTOFDetchid[j]);
            TOTIDIndex.push_back(FTOFDetchid[j]);
            ho[FTOFDetchid[j]] = (TH1D *)hLowU->Clone(buff);
            //
            // *** High TOT distribution ***
            // *********
            c = cdC(6);
            hHighU->GetXaxis()->SetRangeUser(UL, UR);
            hHighU->Draw();

            DrawMyHist(hHighU, "", "", 1, 2, 1);
            DrawMyHist(hHighU_Pick, "", "", 2, 2, 7);
            hHighU_Pick->Draw("SAME");
            leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
            leg->AddEntry(hHighU, "No cut", "l");
            leg->AddEntry(hHighU_Pick, "With cut", "l");
            leg->Draw();

            ll->Draw();
            //sprintf(buff, "%sHighTOTL.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_HighTOTL", cName.Data()));

            //hHighU->GetXaxis()->SetRangeUser(0, 1.5e3);
            //sprintf(buff, "%sHighTOTS.png", PngName.Data());
            //c->SaveAs(buff);

            sprintf(buff, "hho%d", FTOFDetchid[j]);
            hho[FTOFDetchid[j]] = (TH1D *)hHighU->Clone(buff);
            c->Close();

            //
            // *** LowTOT distribution - Risetime ***
            // *********
            c = cdC(7);
            hLowRT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hLowRT, "", "");

            ll->Draw();
            //sprintf(buff, "%sLowRT.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_LowRT", cName.Data()));
            c->Close();

            //
            // *** HighTOT distribution - Risetime ***
            // *********
            c = cdC(8);
            hHighRT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hHighRT, "", "");

            ll->Draw();
            //sprintf(buff, "%sHighRT.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_HighRT", cName.Data()));
            c->Close();

            //
            // *** LowThTime-LowTOT distribution ***
            // *********
            c = cdC(9);
            hLowAT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);
            DrawMy2dHist(hLowAT, "", "");

            ll->Draw();
            //sprintf(buff, "%sLowAT.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_LowAT", cName.Data()));
            c->Close();

            //
            // *** LowThTime-HighTOT distribution ***
            // *********
            c = cdC(10);
            hHighAT->Draw("colz");
            SetMyPad(gPad, 0.15, 0.1, 0.1, 0.15);

            DrawMy2dHist(hHighAT, "", "");

            ll->Draw();
            //sprintf(buff, "%sHighAT.png", PngName.Data());
            //c->SaveAs(buff);
            if (!existflag)
                opfile->WriteObject(c, Form("%s_HighAT", cName.Data()));
            c->Close();

            cout << "loop state:" << LoopOn << ", No. " << j << endl;
            if (!LoopOn)
                break;
        }
    }

    ValidFTOF = eventcounter;
    TString PngName = fileName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5);
    /*

    //
    // ---------draw tot of PMT--------//
    //
    TCanvas *c_all1;
    TCanvas *c_all2;
    TCanvas *c_all3;
    TCanvas *c_all4;
    TCanvas *p1;
    TCanvas *p2;
    TCanvas *p3;
    TCanvas *p4;

    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);

    c_all1->Divide(8, 4);
    c_all2->Divide(8, 4);
    c_all3->Divide(8, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowTOTL", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighTOTL", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_2dU", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowTOT_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighTOT_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%s2dTOT_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);

    //
    // ---------draw time of PMT--------//
    //
    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all1->Divide(8, 4);
    c_all2->Divide(8, 4);
    c_all3->Divide(8, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowThTime", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighThTime", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_RiseTime", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowThtime_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighThtime_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%sRiseTime_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);

    //
    // ---------draw 2d relationship  of PMT--------//
    //
    c_all1 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all2 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all3 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all4 = new TCanvas(Form("c%d", counter++), Form("c%d", counter++), 1600, 1400);
    c_all1->Divide(8, 4);
    c_all2->Divide(8, 4);
    c_all3->Divide(8, 4);
    c_all4->Divide(8, 4);

    for (int i = 0; i < FTOFChNum; i++)
    {
        if (Eabled[i])
        {
            p1 = (TCanvas *)opfile->Get(Form("CH%d_LowRT", FTOFDetchid[i]));
            p2 = (TCanvas *)opfile->Get(Form("CH%d_HighRT", FTOFDetchid[i]));
            p3 = (TCanvas *)opfile->Get(Form("CH%d_LowAT", FTOFDetchid[i]));
            p4 = (TCanvas *)opfile->Get(Form("CH%d_HighAT", FTOFDetchid[i]));

            c_all1->cd(FTOFDetchid[i] + 1);
            p1->DrawClonePad();
            c_all2->cd(FTOFDetchid[i] + 1);
            p2->DrawClonePad();
            c_all3->cd(FTOFDetchid[i] + 1);
            p3->DrawClonePad();
            c_all4->cd(FTOFDetchid[i] + 1);
            p4->DrawClonePad();
        }
    }
    sprintf(buff, "%sLowRT_of_FTOF.png", PngName.Data());
    c_all1->SaveAs(buff);
    sprintf(buff, "%sHighRT_of_FTOF.png", PngName.Data());
    c_all2->SaveAs(buff);
    sprintf(buff, "%sLowAT_of_FTOF.png", PngName.Data());
    c_all3->SaveAs(buff);
    sprintf(buff, "%sHighAT_of_FTOF.png", PngName.Data());
    c_all4->SaveAs(buff);

    //
    // ---------draw NPMT vs TR --------//
    //
    c = cdC(counter++);
    DrawMy2dHist(h2dPMT, "", "", 1, 2);
    //ht->Rebin(2);
    h2dPMT->Draw("colz");
    h2dPMT->GetXaxis()->SetNdivisions(505);
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    //sprintf(buff, "%s2dPMT.png", PngName.Data());
    //c->SaveAs(buff);
    opfile->WriteObject(c, Form("Eff_of_FTOF_vs_Time"));
    //return;

*/

    //
    // ---------draw efficiency of PMT--------//
    //
    counter = 11;

    double PMTEff[FTOFChNum];
    //double PMTEff_Pick[FTOFChNum];
    c = cdC(counter++);
    DrawMyHist(hPMTID, "Channel ID", "Counts", 1, 3);
    //ht->Rebin(2);
    hPMTID->Draw();
    hPMTID->GetXaxis()->SetNdivisions(505);
    c->Update();
    SetEffstats(c, hPMTID, NumofTriggers, RecordTriggers, ValidFTOF);

    //DrawMyHist(hPMTID_Pick, "", "", 2, 2, 7);
    //hPMTID_Pick->Draw("same");
    //leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    //leg->AddEntry(hPMTID, "No cut", "l");
    //leg->AddEntry(hPMTID_Pick, "With cut", "l");
    //leg->Draw();
    /*
    for (int i = 0; i < FTOFChNum; i++)
    {
        if (i > 8)
        {
            ll = DrawMyLatex("... ...", 0.35, 0.7 - 0.05 * i, 42, 0.035, 2);
            ll->Draw("same");
            break;
        }
        PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / RecordTriggers;
        //PMTEff_Pick[i] = hPMTID_Pick->GetBinContent(hPMTID_Pick->FindBin(i)) / RecordTriggers;
        sprintf(buff, "ChID=%d,Ratio=%.1f%%", i, PMTEff[i] * 100);
        //sprintf(buff, "PMTID=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i + 1, PMTEff[i] * 100, PMTEff_Pick[i] * 100);
        ll = DrawMyLatex(buff, 0.3, 0.7 - 0.05 * i, 42, 0.04, 2);
        ll->Draw("same");
    }
    */
    double PMTEffsum = 0;
    int PMTEffctr = 0;
    for (int i = 0; i < FTOFChNum; i++)
    {

        PMTEff[i] = hPMTID->GetBinContent(hPMTID->FindBin(i)) / RecordTriggers;
        if (PMTEff[i] > 0)
        {
            PMTEffsum += PMTEff[i];
            PMTEffctr++;
        }
        //PMTEff_Pick[i] = hPMTID_Pick->GetBinContent(hPMTID_Pick->FindBin(i)) / RecordTriggers;
    }

    sprintf(buff, "Average Ratio=%.1f%%", PMTEffsum / PMTEffctr * 100);
    //sprintf(buff, "PMTID=%d,Eff=%.1f%%,#color[2]{Eff=%.1f%%}", i + 1, PMTEff[i] * 100, PMTEff_Pick[i] * 100);
    ll = DrawMyLatex(buff, 0.3, 0.3, 42, 0.06, 2);
    ll->Draw("same");
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_FTOFCH.png", PngName.Data());
    c->SaveAs(buff);
    if (!existflag)
        opfile->WriteObject(c, Form("Eff_of_FTOFCH"));
    //return;
    //
    // ---------draw hit pos of PMT--------//
    //

    c = cdC(counter++, 1300, 600);
    DrawMy2dHist(hpos, "", "");
    //ht->Rebin(2);
    hpos->Draw("colz");

    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sPos_of_FTOFCH.png", PngName.Data());
    c->SaveAs(buff);
    if (!existflag)
        opfile->WriteObject(hpos, Form("Pos_of_FTOFCH"));

    //return;

    //
    // *** Number of fired PMT ***
    // **
    double Eff[FTOFChNum + 1];
    double Eff_Pick[FTOFChNum + 1];
    TLatex *la;

    c = cdC(counter++);
    //ht->Rebin(2);
    //DrawMyHist(hNPMT_Pick, "", "", 2, 2, 7);
    //hNPMT_Pick->Draw("SAME");
    //leg = DrawMyLeg(0.6, 0.8, 0.9, 0.9);
    //leg->AddEntry(hNPMT, "No cut", "l");
    //leg->AddEntry(hNPMT_Pick, "With cut", "l");
    //leg->Draw();
    int noevent = hNPMT->GetBinContent(hNPMT->FindBin(0));
    if (noevent < 10)
    {
        noevent = RecordTriggers - ValidFTOF;
        hNPMT->SetBinContent(hNPMT->FindBin(0), noevent);
    }
    DrawMyHist(hNPMT, "", "", 1, 3);
    hNPMT->Draw();
    //hNPMT->GetXaxis()->SetNdivisions(505);
    c->Update();
    SetEffstats(c, hNPMT, NumofTriggers, RecordTriggers, ValidFTOF);

    /*
    Eff[0] = hNPMT->GetBinContent(hNPMT->FindBin(0)) / RecordTriggers;
    sprintf(buff, "No ch response,Ratio=%.1f%%", Eff[0] * 100);
    la = DrawMyLatex(buff, 0.25, 0.7, 42, 0.05, 1);
    la->Draw("same");
    //    Eff_Pick[i] = hNPMT_Pick->GetBinContent(hNPMT_Pick->FindBin(i)) / RecordTriggers;

    for (int i = 1; i < FTOFChNum + 1; i++)
    {
        if (i > 8)
        {
            la = DrawMyLatex("... ...", 0.3, 0.7 - 0.05 * i, 42, 0.035);
            la->Draw("same");
            break;
        }

        Eff[i] = hNPMT->GetBinContent(hNPMT->FindBin(i)) / RecordTriggers;
        //Eff_Pick[i] = hNPMT_Pick->GetBinContent(hNPMT_Pick->FindBin(i)) / eventcounter;

        sprintf(buff, "Nch=%d,Ratio=%.1f%%,#color[2]{(Valid)%.1f%%}", i, Eff[i] * 100, Eff[i] * 100 * RecordTriggers / ValidFTOF);
        la = DrawMyLatex(buff, 0.25, 0.7 - 0.05 * i, 42, 0.04);
        la->Draw("same");
    }
    */
    sprintf(buff, "Fired channels=%.1f", hNPMT->GetMean());
    la = DrawMyLatex(buff, 0.4, 0.4, 42, 0.06);
    la->Draw("same");
    //fit = gausfit(ht, 20e-3, 3, 3, 1, tL, tR);
    //sprintf(buff, "TR=%.0fps", fit->GetParameter(2) * 1e3);
    //la = DrawMyLatex(buff, 0.2, 0.4);
    sprintf(buff, "%sEff_of_FTOF.png", PngName.Data());
    c->SaveAs(buff);
    if (!existflag)
        opfile->WriteObject(c, Form("Eff_of_FTOF"));

    //
    // ---------draw LowTOT together --------//
    //
    co = cdC(counter++);
    int max = 0;
    int temp = 0;
    sort(TOTIDIndex.begin(), TOTIDIndex.end());
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {
        ho[TOTIDIndex.at(i)]->Rebin(10);
        cout << "lowTOThist: " << i << endl;
        temp = ho[TOTIDIndex.at(i)]->GetMaximum();
        max = max < temp ? temp : max;
        FillHistArea(ho[TOTIDIndex.at(i)], clr[i % 8], 0.8, fst[2]);
        sprintf(buff, "Ch%d", TOTIDIndex.at(i));
        lego->AddEntry(ho[TOTIDIndex.at(i)], buff, "lpf");
        if (i == 0)
            ho[TOTIDIndex.at(i)]->Draw();

        else
            ho[TOTIDIndex.at(i)]->Draw("same");
    }
    ho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max * 1.1);
    ho[TOTIDIndex.at(0)]->GetXaxis()->SetRangeUser(0, 2.5e3);
    ho[TOTIDIndex.at(0)]->SetStats(0);
    co->Modified();
    lego->Draw();
    sprintf(buff, "%sLowTOTDrawtogether_of_FTOF.png", PngName.Data());
    co->SaveAs(buff);
    if (!existflag)
        opfile->WriteObject(c, Form("all_ch_lowtot"));

    //
    // ---------draw HighTOT together --------//
    //
    lego->Clear();
    co = cdC(counter++);
    max = 0;
    temp = 0;
    for (int i = 0; i < TOTIDIndex.size(); i++)
    {
        hho[TOTIDIndex.at(i)]->Rebin(10);

        cout << "highTOThist: " << i << endl;
        temp = hho[TOTIDIndex.at(i)]->GetMaximum();
        max = max < temp ? temp : max;
        FillHistArea(hho[TOTIDIndex.at(i)], clr[i % 8], 0.1, fst[2]);
        sprintf(buff, "Ch%d", TOTIDIndex.at(i));
        lego->AddEntry(hho[TOTIDIndex.at(i)], buff, "lpf");
        if (i == 0)
            hho[TOTIDIndex.at(i)]->Draw();
        else
            hho[TOTIDIndex.at(i)]->Draw("same");
    }
    hho[TOTIDIndex.at(0)]->GetYaxis()->SetRangeUser(0, max * 1.1);
    hho[TOTIDIndex.at(0)]->GetXaxis()->SetRangeUser(0, 2.5e3);
    ho[TOTIDIndex.at(0)]->SetStats(0);
    co->Modified();
    lego->Draw();
    sprintf(buff, "%sHighTOTDrawtogether__of_FTOF.png", PngName.Data());
    co->SaveAs(buff);
    if (!existflag)
        opfile->WriteObject(c, Form("all_ch_hightot"));
}

void ReadAllData(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210208", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");
    if (fileDir.EndsWith("/"))
        fileDir.Remove(fileDir.Length() - 1, 1);

    TString fileDir1(fileDir + "/Tracker/");
    TString fileDir2(fileDir + "/T0/");
    TString fileDir3(fileDir + "/Trigger/");
    TString fileDir4(fileDir + "/Combine/");
    vector<TString> datList1;
    vector<TString> datList2;
    vector<TString> datList3;

    GetFileList(fileDir1, ".dat", datList1); //for tracker
    GetFileList(fileDir2, ".txt", datList2); //for FTOF
    GetFileList(fileDir3, ".txt", datList3); // for trigger

    if (datList1.size() > 0)
    {
        ReadTrackAGTData2Root(datList1, GenPath(Tracker, RAW, fileDir), force);
        AnalysisAGET(GenPath(Tracker, RAW, fileDir), GenPath(Tracker, PRO, fileDir), force);
    }

    if (datList2.size() > 0)
        ReadFTOFData2Root(datList2, GenPath(T0, RAW, fileDir), force);
    if (datList3.size() > 0)
        ReadTrigger2Root(datList3, GenPath(Trigger, RAW, fileDir), force);

    AnalysisFTOF(GenPath(T0, RAW, fileDir), GenPath(Trigger, RAW, fileDir), GenPath(T0, PRO, fileDir), force);
    FilterFTOF(GenPath(T0, PRO, fileDir), GenPath(T0, COM, fileDir), force);
    //CheckTrackerID(GenPath(Tracker, PRO, fileDir));
    Matchall(GenPath(T0, COM, fileDir), GenPath(Tracker, PRO, fileDir), GenPath(ALL, COM, fileDir), 1);
    //PrintT0(GenPath(ALL, COM, fileDir), GenPath(ALL, PIC, fileDir), 1);
}
void ReadTracker(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210207", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");

    cout << " --> The Input fileName is:" << fileDir << endl;
    //  **
    //  * Read Tracker AGET Data *
    //
    TString Trackerrawname;
    TString Trackerananame;
    TString Trackerdecname;
    TString Trackermatname;
    vector<TString> datList;
    TString fileDir1;
    if (fileDir.EndsWith(".dat"))
    {
        fileDir1.Replace(0, 1000, fileDir, fileDir.Last('.'));
        datList.push_back(fileDir);
        cout << " --> The prename is:" << fileDir1 << endl;
        Trackerrawname = fileDir1 + "Tracker-raw.root";
        Trackerananame = fileDir1 + "Tracker-pro.root";
        Trackerdecname = fileDir1 + "Tracker-dec.root";
        Trackermatname = fileDir1 + "Tracker-mat.root";
    }
    else
    {
        if (fileDir.EndsWith("/"))
            fileDir.Remove(fileDir.Length() - 1, 1);
        fileDir1 = fileDir + "/Tracker/";
        GetFileList(fileDir1, ".dat", datList);
        /*
    ReadTrackAGTData2Root(datList, rootrawname, 0);
    rootpedname = path + "/" + date + "-PED.root";
    GenerateAGETPed(rootrawname, rootpedname);
    ReadPed(rootpedname);
    CalculateEff(rootrawname);
    //AnalysisAGET(rootrawname, rootananame);
    //CalculateEff2(rootananame, roothitname);
    */
        Trackerrawname = fileDir + "/Combine/Tracker-raw.root";
        Trackerananame = fileDir + "/Combine/Tracker-pro.root";
        Trackerdecname = fileDir + "/Combine/Tracker-dec.root";
        Trackermatname = fileDir + "/Combine/Tracker-mat.root";
    }
    // #produce raw-root
    ReadTrackAGTData2Root(datList, Trackerrawname, force);

    // #calculate q,pedmean,pedrms,
    AnalysisAGET(Trackerrawname, Trackerananame, force);

    // #decode AGET chn to MM chn
    //DecodeAGET(Trackerananame, Trackerdecname, 1);
    //
    //Match2("", Trackerdecname, Trackermatname, 1);
    //TrackRebuild(Trackermatname, "01");

    // #Check ID
    CheckTrackerID(Trackerananame);

    // #print tracker information
    //PrintTracker(Trackerananame, 0);
    TrackerEff(Trackerananame);

    bench.Show("full");
};
void AnalysisTracker(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210207", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");

    cout << " --> The Input fileName is:" << fileDir << endl;
    //  **
    //  * Read Tracker AGET Data *
    //
    TString Trackerrawname;
    TString Trackerananame;
    TString Trackerdecname;
    TString Trackermatname;
    vector<TString> datList;
    TString fileDir1;
    if (fileDir.EndsWith(".dat"))
    {
        fileDir1.Replace(0, 1000, fileDir, fileDir.Last('.'));
        datList.push_back(fileDir);
        cout << " --> The prename is:" << fileDir1 << endl;
        Trackerrawname = fileDir1 + "Tracker-raw.root";
        Trackerananame = fileDir1 + "Tracker-pro.root";
        Trackerdecname = fileDir1 + "Tracker-dec.root";
        Trackermatname = fileDir1 + "Tracker-mat.root";
    }
    else
    {
        if (fileDir.EndsWith("/"))
            fileDir.Remove(fileDir.Length() - 1, 1);
        fileDir1 = fileDir + "/Tracker/";
        GetFileList(fileDir1, ".dat", datList);
        /*
    ReadTrackAGTData2Root(datList, rootrawname, 0);
    rootpedname = path + "/" + date + "-PED.root";
    GenerateAGETPed(rootrawname, rootpedname);
    ReadPed(rootpedname);
    CalculateEff(rootrawname);
    //AnalysisAGET(rootrawname, rootananame);
    //CalculateEff2(rootananame, roothitname);
    */
        Trackerrawname = fileDir + "/Combine/Tracker-raw.root";
        Trackerananame = fileDir + "/Combine/Tracker-pro.root";
        Trackerdecname = fileDir + "/Combine/Tracker-dec.root";
        Trackermatname = fileDir + "/Combine/Tracker-mat.root";
    }
    // #produce raw-root
    ReadTrackAGTData2Root(datList, Trackerrawname, force);

    // #calculate q,pedmean,pedrms,
    AnalysisAGET(Trackerrawname, Trackerananame, force);

    bench.Show("full");
};
void GenTrackerdecdata(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210207", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");

    cout << " --> The Input fileName is:" << fileDir << endl;
    //  **
    //  * Read Tracker AGET Data *
    //
    TString Trackerrawname;
    TString Trackerananame;
    TString Trackerdecname;
    TString Trackermatname;
    vector<TString> datList;
    TString fileDir1;
    if (fileDir.EndsWith(".dat"))
    {
        fileDir1.Replace(0, 1000, fileDir, fileDir.Last('.'));
        datList.push_back(fileDir);
        cout << " --> The prename is:" << fileDir1 << endl;
        Trackerrawname = fileDir1 + "Tracker-raw.root";
        Trackerananame = fileDir1 + "Tracker-pro.root";
        Trackerdecname = fileDir1 + "Tracker-dec.root";
        Trackermatname = fileDir1 + "Tracker-mat.root";
    }
    else
    {
        if (fileDir.EndsWith("/"))
            fileDir.Remove(fileDir.Length() - 1, 1);
        fileDir1 = fileDir + "/Tracker/";
        GetFileList(fileDir1, ".dat", datList);
        /*
    ReadTrackAGTData2Root(datList, rootrawname, 0);
    rootpedname = path + "/" + date + "-PED.root";
    GenerateAGETPed(rootrawname, rootpedname);
    ReadPed(rootpedname);
    CalculateEff(rootrawname);
    //AnalysisAGET(rootrawname, rootananame);
    //CalculateEff2(rootananame, roothitname);
    */
        Trackerrawname = fileDir + "/Combine/Tracker-raw.root";
        Trackerananame = fileDir + "/Combine/Tracker-pro.root";
        Trackerdecname = fileDir + "/Combine/Tracker-dec.root";
        Trackermatname = fileDir + "/Combine/Tracker-mat.root";
    }
    // #produce raw-root
    ReadTrackAGTData2Root(datList, Trackerrawname, force);

    // #calculate q,pedmean,pedrms,
    AnalysisAGET(Trackerrawname, Trackerananame, force);

    // #decode AGET chn to MM chn
    DecodeAGET4Align(Trackerananame, Trackerdecname, 1);

    bench.Show("full");
};
void ReadT0(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210208", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");
    if (fileDir.EndsWith("/"))
        fileDir.Remove(fileDir.Length() - 1, 1);

    TString fileDir2(fileDir + "/T0/");
    TString fileDir3(fileDir + "/Trigger/");
    TString fileDir4(fileDir + "/Combine/");
    vector<TString> datList2;
    vector<TString> datList3;

    GetFileList(fileDir2, ".txt", datList2); //for FTOF
    GetFileList(fileDir3, ".txt", datList3); // for trigger

    // #produce raw-root
    if (datList2.size() > 0)
        ReadFTOFData2Root(datList2, GenPath(T0, RAW, fileDir), force);
    if (datList3.size() > 0)
        ReadTrigger2Root(datList3, GenPath(Trigger, RAW, fileDir), force);

    // #combine T0 with Trigger
    AnalysisFTOF(GenPath(T0, RAW, fileDir), GenPath(Trigger, RAW, fileDir), GenPath(T0, PRO, fileDir), force);

    // #Pickup valid data and Change data format : combine each channel data to a event
    FilterFTOF(GenPath(T0, PRO, fileDir), GenPath(T0, COM, fileDir), force);

    // #print T0 information
    //PrintT0(GenPath(T0, COM, fileDir), GenPath(T0, PIC, fileDir), force);
    GetTriggerNum(GenPath(Trigger, RAW, fileDir));
    PrintT0(GenPath(T0, COM, fileDir), GenPath(T0, PIC, fileDir), force);
}
void ReadFTOF(TString fileDir = "/mnt/d/Experiment/DATA/cosmicray/20210208", int force = 0)
{
    TBenchmark bench;
    bench.Start("full");
    if (fileDir.EndsWith("/"))
        fileDir.Remove(fileDir.Length() - 1, 1);

    TString fileDir2(fileDir + "/T0/");
    TString fileDir3(fileDir + "/Trigger/");
    TString fileDir4(fileDir + "/Combine/");
    vector<TString> datList2;
    vector<TString> datList3;

    GetFileList(fileDir2, ".txt", datList2); //for FTOF
    GetFileList(fileDir3, ".txt", datList3); // for trigger

    // #produce raw-root
    if (datList2.size() > 0)
        ReadFTOFData2Root(datList2, GenPath(FTOF, RAW, fileDir), force);
    if (datList3.size() > 0)
        ReadTrigger2Root(datList3, GenPath(Trigger, RAW, fileDir), force);

    // #combine FTOF with Trigger
    AnalysisFTOF(GenPath(FTOF, RAW, fileDir), GenPath(Trigger, RAW, fileDir), GenPath(FTOF, PRO, fileDir), force);

    // #Pickup valid data and Change data format : combine each channel data to a event
    FilterFTOF(GenPath(FTOF, PRO, fileDir), GenPath(FTOF, COM, fileDir), force);

    // #print T0 information
    //PrintT0(GenPath(T0, COM, fileDir), GenPath(T0, PIC, fileDir), force);
    GetTriggerNum(GenPath(Trigger, RAW, fileDir));
    PrintFTOF(GenPath(FTOF, COM, fileDir), GenPath(FTOF, PIC, fileDir), force);
}

#if 0
void drawdecoderesults(TString fDecName)
{
    //* Tracker data
    MMHit *MMhit = new MMHit;
    TFile *trackerfile = new TFile(fDecName, "read");
    TTree *trackertree[4];
    int allmin = 9999999;
    int allmax = 0;
    for (int s = 0; s < 4; s++)
    {
        trackertree[s] = (TTree *)trackerfile->Get(Form("tree%d", s));
        trackertree[s]->SetBranchAddress("MMhit", &MMhit);
        cout << "Entry Detector " << s << " ..." << endl;
        int nEntries = trackertree[s]->GetEntriesFast();
    }
}
void drawhitmap(TString fMatName)
{
    CRHit *fhit = new CRHit;
    TFile *file;
    TTree *tree;
    file = new TFile(fMatName, "read");
    tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("CRhit", &fhit);
    int nEntries = tree->GetEntries();
    cout << fMatName << " has Entries: " << nEntries << endl;
    TH2F *fhitmapReal = new TH2F(Form("fhitmapReal"), Form("Real Position hit map for %s", detName[0]), 800, -800 / 2 * 0.4, 800 / 2 * 0.4, 800, -800 / 2 * 0.4, 800 / 2 * 0.4);
    int nbin = 400;
    TH2F *hitmap = new TH2F(Form("fhitmap"), Form("hit map for %s", detName[0]), nbin, 0, nbin, nbin, 0, nbin);
    for (int i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);
        int nMMs = fhit->Tracker.size();
        for (int s = 0; s < nMMs; s++)
        {
            if (fhit->Tracker[s].MMID == 0)
            {
                double xqmax = -1;
                double xmmax = -1;
                for (int ii = 0; ii < (int)fhit->Tracker[ii].MMchn[X].size(); ii++)
                {
                    xmmax = (xqmax < fhit->Tracker[s].MMchn[X][ii].q) ? fhit->Tracker[s].MMchn[X][ii].clustermean[0] : xmmax;
                    xqmax = (xqmax < fhit->Tracker[s].MMchn[X][ii].q) ? fhit->Tracker[s].MMchn[X][ii].q : xqmax;
                }
                double yqmax = -1;
                double ymmax = -1;
                for (int jj = 0; jj < (int)fhit->Tracker[s].MMchn[Y].size(); jj++)
                {
                    ymmax = (yqmax < fhit->Tracker[s].MMchn[Y][jj].q) ? fhit->Tracker[s].MMchn[Y][jj].clustermean[1] : ymmax;
                    yqmax = (yqmax < fhit->Tracker[s].MMchn[Y][jj].q) ? fhit->Tracker[s].MMchn[Y][jj].q : yqmax;
                }
                if (xmmax != -1 && yqmax != -1)
                    hitmap->Fill(xmmax, ymmax);
            }
        }
    }
    hitmap->Draw();
    hitmap->SetMarkerColor(2);
    hitmap->SetMarkerSize(1);
    hitmap->SetMarkerStyle(24);
}

bool TrackRebuild(TString fMatName, TString Runnum)
{
    TString fpath = GetPath(fMatName);
    fName = "RUN" + Runnum;
    ReadOffset();
    Definehist();
    CRHit *fhit = new CRHit;
    TFile *file;
    TTree *tree;
    file = new TFile(fMatName, "read");
    tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("CRhit", &fhit);
    int nEntries = tree->GetEntries();
    cout << fMatName << " has Entries: " << nEntries << endl;

    for (int i = 0; i < nEntries; i++)
    {
        tree->GetEntry(i);
        for (int ii = 0; ii < NumMM; ii++)
        {
            FillTrackHistogram(fhit->Tracker, ii);
        }
        // cout <<  " Trackhist has been filled " << endl;

        for (int xy = 0; xy < 2; xy++)
        {

            //AnalysisTrack(vector<int>{MM0id, MM1id}, MM3id, fhit->Tracker, xy);
            //AnalysisTrack(vector<int>{MM0id, MM3id}, MM1id, fhit->Tracker, xy);
            //AnalysisTrack(vector<int>{MM1id, MM3id}, MM0id, fhit->Tracker, xy);
            AnalysisTrack(vector<int>{MM0id, MM1id, MM2id}, MM3id, fhit->Tracker, xy);
            AnalysisTrack(vector<int>{MM0id, MM1id, MM3id}, MM2id, fhit->Tracker, xy);
            AnalysisTrack(vector<int>{MM0id, MM3id, MM2id}, MM1id, fhit->Tracker, xy);
            AnalysisTrack(vector<int>{MM3id, MM1id, MM2id}, MM0id, fhit->Tracker, xy);
        }
        //   cout << " Track has been analyzed " << endl;
    }

    //--------------------------
    // 画图
    gStyle->SetOptFit(1);

    for (int i = 0; i < NumMM; i++)
    {

        DrawTracker(fpath, i);
    }
    DrawDeviation("c5", vector<int>{MM0id, MM1id, MM2id, MM3id});
    PrintInfo();
    return true;
}
#endif

/*

void Combineall(TString fAnaName1, TString fAnaName2, TString fComName, int force)
{
    if (force != 1)
    {
        FileStat_t fStat;
        gSystem->GetPathInfo(fComName, fStat);
        if (fStat.fSize != 0)
        {
            cout << "ERROR! " << fComName << " is exist!" << endl;
            return;
        }
    }

    TFile *comFile = new TFile(fComName, "recreate");
    if (!comFile->IsOpen())
    {
        cout << "ERROR! " << fComName << " cant open!" << endl;
        return;
    }

    TimeData TimeData;
    AGETHit trackerhit;
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("TimeData", &TimeData);

    //* FTOF data
    TFile *ftoftreefile = new TFile(fAnaName1, "read");
    TTree *ftoftree = (TTree *)ftoftreefile->Get("htree");
    ftoftree->SetBranchAddress("dethit", &TimeData);

    //* Tracker data
    TFile *trackerfile = new TFile(fAnaName2, "read");
    TTree *trackertree[4];

    vector<int> ftofidvec;
    vector<TimeData> ftofhitlist;
    int N;
    N = ftoftree->GetEntriesFast();
    cout << " FTOF entries: " << N << endl;
    for (int i = 0; i < N; i++)
    {
        ftoftree->GetEntry(i);
        ftofhitlist.push_back(TimeData);
        ftofidvec.push_back(TimeData.event);
    }
    int idstart = TMath::MinElement(ftofidvec.size(), &ftofidvec.at(0));
    int idend = TMath::MaxElement(ftofidvec.size(), &ftofidvec.at(0));

    vector<int> trackeridvec[4];
    vector<AGETHit> trackerhitlist[4];
    for (int s = 0; s < 4; s++)
    {
        trackertree[s] = (TTree *)trackerfile->Get(Form("htree%d", s));
        trackertree[s]->SetBranchAddress("dethit", &trackerhit);
        N = trackertree[s]->GetEntriesFast();
        cout << " Tracker entries: " << N << endl;
        for (int i = 0; i < N; i++)
        {
            trackertree[s]->GetEntry(i);
            trackerhitlist[s].push_back(trackerhit);
            trackeridvec[s].push_back(trackerhit.event);
        }
        idstart = idstart < TMath::MinElement(trackeridvec[s].size(), &trackeridvec[s].at(0)) ? idstart : TMath::MinElement(trackeridvec[s].size(), &trackeridvec[s].at(0));
        idend = idend < TMath::MaxElement(trackeridvec[s].size(), &trackeridvec[s].at(0)) ? idend : TMath::MaxElement(trackeridvec[s].size(), &trackeridvec[s].at(0));
    }
    int ftofcounter = 0;
    int trackercounter = 0;

    vector<int> ftofidindexvec;
    vector<int> trackeridindexvec[4];
    int pos = 0;
    TimeData temp1;
    AGETHit temp2;
    vector<TimeData> ftofbuffer;
    vector<AGETHit> trackerbuffer[4];
    for (int i = idstart; i < idend + 1; i++)
    {
        ftofcounter = SearchIDIndex(i, ftofidvec, ftofidindexvec);

        for (int j = 0; j < 4; j++)
        {
            if (SearchIDIndex(i, trackeridvec[j], trackeridindexvec[j]))
                trackercounter++;
        }

        if (ftofcounter == 1 && trackercounter >= 2)
        {

            for (int k = 0; k < ftofidindexvec.size(); k++)
            {
                pos = ftofidindexvec.at(k);
                temp1 = ftofhitlist.at(pos);
                if (temp1.triggertimediff > -700e3 && temp1.triggertimediff < -640e3)
                {

                    ftofbuffer.push_back(temp1);
                }
            }
            checkrepetition(ftofbuffer);

            for (int j = 0; j < 4; j++)
            {
                for (int k = 0; k < trackeridindexvec[j].size(); k++)
                {
                    pos = trackeridindexvec[j].at(k);
                    temp2 = trackerhitlist[j].at(pos);
                    if (temp2.wave_max > 7 * temp2.ped_rms && temp2.q > 5)
                    {
                        trackerbuffer[j].push_back(temp2);
                    }
                }
            }
            trackercounter = CalculatePos(trackerbuffer);
        }
        if (ftofbuffer.size() > 0 && trackercounter > 1)
        {
            for (int i = 0; i < ftofbuffer.size(); i++)
            {
                TimeData = ftofbuffer.at(i);
                tree->Fill();
            }
        }
    }
    comFile->WriteTObject(tree);
    comFile->Flush();
    comFile->Close();
    ftoftreefile->Close();
    trackerfile->Close();
}

void AnalysisAGET(TString fRawName, TString fAnaName)
{
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;

    TFile *rawFile = new TFile(fRawName);
    TFile *anaFile = new TFile(fAnaName, "recreate");
    TTree *fTree = new TTree("ftree", "ftree");
    double wave_max;
    double charge;
    double ped;
    double ped_rms;
    double factor = 2e-7 / 50 / 1e3 / 1e-9; // factor that converts to nC.
    fTree->Branch("event", &event);
    fTree->Branch("board", &board);
    fTree->Branch("chip", &chip);
    fTree->Branch("channel", &channel);
    fTree->Branch("wave_max", &wave_max);
    fTree->Branch("charge", &charge);
    fTree->Branch("ped", &ped);
    fTree->Branch("ped_rms", &ped_rms);
    cout << "--> Analyzing waveform files from " << fRawName << endl;
    TH1D *hNoise = new TH1D("hNoise", "hNoise", 1000, 0, 1000);
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);
    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);
        hNoise->Reset();
        charge = 0;
        ped = 0;
        ped_rms = 0;
        wave_max = 0;
        for (int j = 0; j < (int)wave->size(); j++)
        {
            if ((j > 5 && j < 150) || (j > 320 && j < 500))
            {

                hNoise->Fill(wave->at(j));
            }
            if (j >= 150 && j < 320)
            {
                charge += wave->at(j) * factor;
            }
        }
        ped = hNoise->GetMean();
        ped_rms = hNoise->GetRMS();
        if (hNoise->GetEntries() > 0)
        {
            if (ped - 3 * ped_rms < 50)
                hNoise->Fit("gaus", "RQ", "", 50, ped + 3 * ped_rms);
            else
                hNoise->Fit("gaus", "RQ", "", ped - 3 * ped_rms, ped + 3 * ped_rms);
            ped = hNoise->GetFunction("gaus")->GetParameter(1);
            ped_rms = hNoise->GetFunction("gaus")->GetParameter(2);
        }
        charge = charge - ped * (320 - 150) * factor;
        wave_max = TMath::MaxElement(320 - 150, &(wave->at(150))) - ped;
        fTree->Fill();
    }
    fTree->Write();
    anaFile->Close();
    rawFile->Close();
    cout << "--> Analyzed results has been stored to " << fAnaName << endl
         << endl;
}

void CalculateEff2(TString fAnaName, TString fHitName)
{

    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    double wave_max;
    double charge;
    double ped;
    double ped_rms;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;

    TBranch *b_wave_max;
    TBranch *b_charge;
    TBranch *b_ped;
    TBranch *b_ped_rms;
    TFile *hitFile = new TFile(fHitName, "recreate");
    TTree *htree[4];
    AGETHit dethit;
    for (int i = 0; i < 4; i++)
    {

        htree[i] = new TTree(Form("htree%d", i), Form("htree%d", i));
        htree[i]->Branch("dethit", &dethit);
    }

    TFile *anaFile = new TFile(fAnaName, "read");
    TTree *tree = (TTree *)anaFile->Get("ftree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave_max", &wave_max, &b_wave_max);
    tree->SetBranchAddress("charge", &charge, &b_charge);
    tree->SetBranchAddress("ped", &ped, &b_ped);
    tree->SetBranchAddress("ped_rms", &ped_rms, &b_ped_rms);
    int eventcounter[8];
    int nentries = tree->GetEntriesFast();

    int boardNo = boardList.size();
    int chipNo = chipList.size();
    int dim = boardNo * chipNo * 64;
    int chNo = 64; // 64*6=384

    vector<TH1I *> hevent;
    hevent.resize(dim);
    TH1I *hid = new TH1I("hid", "hid", 100e4, 1, 100e4 + 1);
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                TH1I *tmp = (TH1I *)gDirectory->Get(Form("Event_%d_%d_%d", boardList[i], chipList[j], k));
                if (tmp != NULL)
                    delete tmp;
                hevent[i * (64 * 4) + j * 64 + k] = new TH1I(Form("Event_%d_%d_%d", boardList[i], chipList[j], k), Form("Event for Board=%d, Chip=%d, Channel=%d", boardList[i], chipList[j], k), 100e4, 1, 100e4 + 1);
            }

    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);
        hid->Fill(event);
        dethit.event = 0;
        dethit.q = 0;
        dethit.wave_max = 0;
        dethit.ped = 0;
        dethit.ped_rms = 0;
        dethit.hit.first = 0;
        dethit.hit.second = 0;
        if (wave_max > 7 * ped_rms && charge > 5)
        {
            hevent[IndexBD(board) * 64 * 4 + IndexChip(chip) * 64 + channel]->Fill(event);
            dethit.event = event;
            dethit.q = charge;
            dethit.wave_max = wave_max;
            dethit.ped = ped;
            dethit.ped_rms = ped_rms;
            if (IndexBD(board) == 0)
            {
                // y axis
                dethit.hit.first = -999;
                dethit.hit.second = channel;
            }
            else
            {
                //x axis
                dethit.hit.first = channel;
                dethit.hit.second = -999;
            }

            htree[IndexChip(chip)]->Fill();
        }
    }
    for (int i = 0; i < 4; i++)
    {
        hitFile->WriteTObject(htree[i]);
    }

    RecordStart = hid->FindFirstBinAbove(1);
    RecordEnd = hid->FindLastBinAbove(1) + 1;

    TH1I *heff[3];
    TH1I *hchipeff = new TH1I("hchipeff", "", 9, 0, 9);
    DrawMyHist(hchipeff, "chipNum", "Counts", 1, 3, 1);
    TH1I *hDeteff[4];
    TH2I *hDetHit[4];
    TH1I *hXYHit[8];
    int channelcounter;
    int boardcounter[2];
    int totalcounter;
    Color_t clr[] = {2, 4, 1};
    Style_t stl[] = {7, 7, 1};
    TString boardtype[] = {"Board1", "Board2", "BoardALL"};
    TLegend *leg = DrawMyLeg(0.7, 0.4, 0.85, 0.55);
    TLatex *la;
    for (int i = 0; i < 3; i++)
    {
        heff[i] = new TH1I(Form("heff%d", i), "", 9, 0, 9);
        DrawMyHist(heff[i], "Nchips", "Counts", clr[i], 3, stl[i]);
        leg->AddEntry(heff[i], Form("%s", boardtype[i].Data()), "lp");
    }
    for (int i = 0; i < 4; i++)
    {
        //x-2,y-1
        hDeteff[i] = new TH1I(Form("heff%d", i), "", 3, 0, 3);
        //leg->AddEntry(hDeteff[i], Form("Det%d",i), "lp");

        hDetHit[i] = new TH2I(Form("hDetHit%d", i), "", chNo, 0, chNo, chNo, 0, chNo);
        hXYHit[4 + i] = new TH1I(Form("HitX%d", i), "X Strip hits for single event", chNo, 0, chNo);
        hXYHit[i] = new TH1I(Form("HitY%d", i), "Y Strip hits for single event", chNo, 0, chNo);
    }
    vector<int> xyhit[8];
    pair<double, double> hit;
    for (int s = RecordStart; s < RecordEnd; s++)
    {

        boardcounter[0] = 0;
        boardcounter[1] = 0;
        totalcounter = 0;
        memset(xyhit, 0, sizeof(xyhit));
        for (int i = 0; i < boardNo; i++)
        {

            for (int j = 0; j < chipNo; j++)

            {
                channelcounter = 0;
                for (int k = 0; k < 64; k++)
                {

                    if (hevent[i * 64 * 4 + j * 64 + k]->GetBinContent(s) > 0)
                    {

                        channelcounter++;
                        hXYHit[i * 4 + j]->Fill(k);
                        xyhit[i * 4 + j].push_back(k);
                    }
                }
                if (channelcounter > 0)
                {
                    hDeteff[j]->Fill(i + 1);
                    hchipeff->Fill(i * 4 + j);
                    boardcounter[i]++;
                    totalcounter++;
                }
            }
        }
        for (int i = 0; i < 4; i++)
        {
            for (int k = 0; k < (int)xyhit[i].size(); k++)
                for (int j = 0; j < (int)xyhit[i + 4].size(); j++)
                {
                    hDetHit[i]->Fill(xyhit[4 + i].at(j), xyhit[i].at(k));
                }
        }
        //cout<<"counter1="<<counter[0]<<", counter2="<<counter[1]<<", counter3="<<counter[2]<<endl;

        heff[0]->Fill(boardcounter[0]);
        heff[1]->Fill(boardcounter[1]);
        heff[2]->Fill(totalcounter);
    }

    TCanvas *c = cdC(0);
    heff[1]->Draw();
    heff[0]->Draw("same");
    heff[2]->Draw("same");
    leg->Draw();

    double Eff[8];
    RecordTriggers = heff[2]->Integral(2, 100);
    cout << "--> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl;
    cout << "RecordTriggers= " << RecordTriggers << endl;
    for (int i = 0; i < 8; i++)
    {

        Eff[i] = heff[2]->GetBinContent(i + 2) / RecordTriggers;
        sprintf(buff, "NChip=%d,Eff=%.1f%%", i + 1, Eff[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    TString PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "NChip.png";
    c->SaveAs(PngName);
    c = cdC(1);
    hchipeff->Draw();
    //leg->Draw();

    double ChipEff[8];
    TString XY[] = {"Y", "X"};
    for (int i = 0; i < 8; i++)
    {

        ChipEff[i] = hchipeff->GetBinContent(i + 1) / RecordTriggers;
        sprintf(buff, "ChipNum=%d,MM%d%s,Eff=%.1f%%", i, i % 4 + 1, XY[i / 4].Data(), ChipEff[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "EveryChip.png";
    c->SaveAs(PngName);

    c = new TCanvas("c3", "c3", 1500, 1200);
    c->Divide(4, 4);
    double detEff[2];
    for (int i = 0; i < 4; i++)
    {
        c->cd(i + 1);
        detEff[0] = 0;
        detEff[1] = 0;
        hDeteff[i]->Draw();
        for (int j = 0; j < 2; j++)
        {

            detEff[j] = hDeteff[i]->GetBinContent(j + 2) / RecordTriggers;
            sprintf(buff, "MM%d%s,Eff=%.1f%%", i, XY[j].Data(), detEff[j] * 100);
            la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * j, 42, 0.05);
            la->Draw("same");
        }
        c->cd(5 + i);
        hXYHit[i]->Draw();
        c->cd(9 + i);
        hXYHit[4 + i]->Draw();
        c->cd(13 + i);
        //hDetHit[i]->SetMarkerStyle(4);
        hDetHit[i]->Draw("colz");
    }
    PngName = fAnaName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "HitDistributed.png";
    c->SaveAs(PngName);
    //hitFile->Close();
    //anaFile->Close();
}

void GenerateAGETPed(TString fRawName, TString fPedName)
{
    //init
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    vector<TH1D *> hNoise;
    TH1I *hid = new TH1I("hid", "hid", 100e4, 1, 100e4 + 1);
    //int boardNo = 3;
    int boardNo = 2;
    int chipNo = 4;
    int dim = boardNo * chipNo * 64;

    hNoise.resize(dim);
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", boardList[i], chipList[j], k));
                if (tmp != NULL)
                    delete tmp;
                hNoise[i * (64 * 4) + j * 64 + k] = new TH1D(Form("Ped_%d_%d_%d", boardList[i], chipList[j], k), Form("Pedestal for Board=%d, Chip=%d, Channel=%d", boardList[i], chipList[j], k), 1000, 0, 1000);
            }

    //read raw file
    TFile *rawFile = new TFile(fRawName);
    TFile *pedFile = new TFile(fPedName, "recreate");

    cout << "--> Analyzing pedestal files from " << fRawName << endl;
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);

    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);
        hid->Fill(event);
        for (int j = 0; j < (int)wave->size(); j++)
        {
            if ((j > 5 && j < 150) || (j > 320 && j < 500))
            {

                hNoise[(board - 1) * (64 * 4) + (chip - 10) * 64 + channel]->Fill(wave->at(j));
            }
        }
    }
    hid->Write();

    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {

                double mean = hNoise[i * (64 * 4) + j * 64 + k]->GetMean();
                double rms = hNoise[i * (64 * 4) + j * 64 + k]->GetRMS();
                int id = i * (64 * 4) + j * 64 + k;
                if (hNoise[id]->GetEntries() > 0)
                {
                    if (mean - 3 * rms < 50)
                        hNoise[id]->Fit("gaus", "RQ", "", 50, mean + 3 * rms);
                    else
                        hNoise[id]->Fit("gaus", "RQ", "", mean - 3 * rms, mean + 3 * rms);
                }
                hNoise[id]->Write();
            }

    pedFile->Write();
    pedFile->Close();
    rawFile->Close();
    cout << "--> Pedestal has been stored to " << fPedName << endl
         << endl;
}
bool ReadPed(TString fPedName)
{

    TFile *fPedFile = new TFile(fPedName);
    if (!fPedFile->IsOpen())
    {
        cout << "#### Can't open pedestal " << fPedName << " to read. Please check the path." << endl;
        return false;
    }
    cout << "--> Now reading the pedestal file: " << fPedName << endl;

    int boardNo = boardList.size();
    int chipNo = chipList.size();
    int dim = boardNo * chipNo * 64;

    SetNPedMean();
    SetNPedRMS();

    for (int ii = 0; ii < boardNo; ii++)
        for (int jj = 0; jj < chipNo; jj++)
            for (int kk = 0; kk < 64; kk++)
            {
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Ped_%d_%d_%d", boardList[ii], chipList[jj], kk));
                if (tmp == 0)
                {
                    cout << "Check the TrackAGET pedestal file: Ped_" << boardList[ii] << "_" << chipList[jj] << "_" << kk << " exists or not" << endl;
                    continue;
                }
                double mean = 0;
                double rms = 0;
                if (tmp->GetEntries() > 0)
                {
                    mean = tmp->GetFunction("gaus")->GetParameter(1);
                    rms = tmp->GetFunction("gaus")->GetParameter(2);
                }
                SetPedMean(ii, jj, kk, mean);
                SetPedRMS(ii, jj, kk, rms);
            }
    TH1I *tmp = (TH1I *)gDirectory->Get("hid");
    RecordStart = tmp->FindFirstBinAbove(1);
    RecordEnd = tmp->FindLastBinAbove(1) + 1;

    cout << "--> Pedestal file has been read.\n"
         << endl;

    return true;
}


void CalculateEff(TString fRawName)
{

    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    vector<int> *pos = 0;
    double wave_max;

    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;
    TBranch *b_pos;
    TBranch *b_wave_max;

    TFile *rawFile = new TFile(fRawName, "read");
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);
    tree->SetBranchAddress("pos", &pos, &b_pos);
    tree->SetBranchAddress("wave_max", &wave_max, &b_wave_max);

    int eventcounter[8];
    int nentries = tree->GetEntriesFast();

    int boardNo = boardList.size();
    int chipNo = chipList.size();
    int dim = boardNo * chipNo * 64;

    vector<TH1I *> hevent;
    hevent.resize(dim);
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)

        {
            TH1I *tmp = (TH1I *)gDirectory->Get(Form("Event_%d_%d", boardList[i], chipList[j]));
            if (tmp != NULL)
                delete tmp;
            hevent[i * 4 + j] = new TH1I(Form("Event_%d_%d", boardList[i], chipList[j]), Form("Event for Board=%d, Chip=%d", boardList[i], chipList[j]), RecordEnd - RecordStart, RecordStart, RecordEnd);
        }

    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);
        if (wave_max - GetPedMean(board, chip, channel) > 7 * GetPedRMS(board, chip, channel))
        {
            hevent[IndexBD(board) * 4 + IndexChip(chip)]->Fill(event);
        }
    }
    TH1I *heff[3];
    TH1I *hchipeff = new TH1I("hchipeff", "", 9, 0, 9);
    DrawMyHist(hchipeff, "chipNum", "Counts", 1, 3, 1);
    int counter[3];
    Color_t clr[] = {2, 4, 1};
    Style_t stl[] = {7, 7, 1};
    TString boardtype[] = {"Board1", "Board2", "BoardALL"};
    TLegend *leg = DrawMyLeg(0.7, 0.4, 0.85, 0.55);
    TLatex *la;
    for (int i = 0; i < 3; i++)
    {
        counter[i] = 0;
        heff[i] = new TH1I(Form("heff%d", i), "", 9, 0, 9);
        DrawMyHist(heff[i], "Nchips", "Counts", clr[i], 3, stl[i]);
        leg->AddEntry(heff[i], Form("%s", boardtype[i].Data()), "lp");
    }

    for (int k = RecordStart; k < RecordEnd; k++)
    {

        counter[0] = 0;
        counter[1] = 0;
        counter[2] = 0;
        for (int i = 0; i < boardNo; i++)
        {

            for (int j = 0; j < chipNo; j++)

            {
                if (hevent[i * 4 + j]->GetBinContent(k) > 0)
                {
                    hchipeff->Fill(i * 4 + j);
                    counter[2]++;
                    counter[i]++;
                }
            }
        }
        //cout<<"counter1="<<counter[0]<<", counter2="<<counter[1]<<", counter3="<<counter[2]<<endl;

        heff[0]->Fill(counter[0]);
        heff[1]->Fill(counter[1]);
        heff[2]->Fill(counter[2]);
    }

    TCanvas *c = cdC(0);
    heff[1]->Draw();
    heff[0]->Draw("same");
    heff[2]->Draw("same");
    leg->Draw();

    double Eff[8];
    RecordTriggers = heff[2]->Integral(2, 100);
    cout << "--> Trigger Start = " << RecordStart << ", Trigger End =" << RecordEnd << endl;
    cout << "RecordTriggers= " << RecordTriggers << endl;
    for (int i = 0; i < 8; i++)
    {

        Eff[i] = heff[2]->GetBinContent(i + 2) / RecordTriggers;
        sprintf(buff, "NChip=%d,Eff=%.1f%%", i + 1, Eff[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    TString PngName = fRawName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "NChip.png";
    c->SaveAs(PngName);
    c = cdC(1);
    hchipeff->Draw();
    //leg->Draw();

    double ChipEff[8];
    TString XY[] = {"Y", "X"};
    for (int i = 0; i < 8; i++)
    {

        ChipEff[i] = hchipeff->GetBinContent(i + 1) / RecordTriggers;
        sprintf(buff, "ChipNum=%d,MM%d%s,Eff=%.1f%%", i, i % 4 + 1, XY[i / 4].Data(), ChipEff[i] * 100);
        la = DrawMyLatex(buff, 0.3, 0.3 + 0.07 * i, 42, 0.05);
        la->Draw("same");
    }
    PngName = fRawName.Copy();
    PngName = PngName.Remove(PngName.Length() - 5, 5) + "EveryChip.png";
    c->SaveAs(PngName);
}

void calculateQ(TString fRawName, TString fQName)
{
    int event;
    UShort_t board;
    UShort_t chip;
    UShort_t channel;
    vector<double> *wave = 0;
    TBranch *b_event;
    TBranch *b_board;
    TBranch *b_chip;
    TBranch *b_channel;
    TBranch *b_wave;

    vector<TH1D *> hQ;
    double charge[2 * 4 * 64];
    double factor = 2e-7 / 50 / 1e3 / 1e-12; // factor that converts to pC.

    int boardNo = 2;
    int chipNo = 4;
    int dim = boardNo * chipNo * 64;
    hQ.resize(dim);
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {
                TH1D *tmp = (TH1D *)gDirectory->Get(Form("Q_%d_%d_%d", boardList[i], chipList[j], k));
                if (tmp != NULL)
                    delete tmp;
                hQ[i * (64 * 4) + j * 64 + k] = new TH1D(Form("Q_%d_%d_%d", boardList[i], chipList[j], k), Form("Charge for Board=%d, Chip=%d, Channel=%d", boardList[i], chipList[j], k), 1000, 0, 1000);
            }

    TFile *rawFile = new TFile(fRawName);
    TFile *QFile = new TFile(fQName, "recreate");

    cout << "--> Analyzing pedestal files from " << fRawName << endl;
    TTree *tree = (TTree *)rawFile->Get("tree");
    tree->SetBranchAddress("event", &event, &b_event);
    tree->SetBranchAddress("board", &board, &b_board);
    tree->SetBranchAddress("chip", &chip, &b_chip);
    tree->SetBranchAddress("channel", &channel, &b_channel);
    tree->SetBranchAddress("wave", &wave, &b_wave);
    int nentries = tree->GetEntriesFast();
    for (int i = 0; i < nentries; i++)
    {
        Long64_t ii = tree->LoadTree(i);
        if (ii < 0)
            break;

        tree->GetEntry(ii);
        for (int j = 150; j < 320; j++)
        {
            charge[IndexBD(board) * (64 * 4) + IndexChip(chip) * 64 + channel] += (wave->at(j) - GetPedMean(board, chip, channel)) * factor;
        }
        hQ[IndexBD(board) * (64 * 4) + IndexChip(chip) * 64 + channel]->Fill(charge[IndexBD(board) * (64 * 4) + IndexChip(chip) * 64 + channel]);
    }
    for (int i = 0; i < boardNo; i++)
        for (int j = 0; j < chipNo; j++)
            for (int k = 0; k < 64; k++)
            {

                hQ[i * (64 * 4) + j * 64 + k]->Write();
            }

    QFile->Write();
    QFile->Close();
    rawFile->Close();
    cout << "--> Pedestal has been stored to " << fQName << endl
         << endl;
}
*/
