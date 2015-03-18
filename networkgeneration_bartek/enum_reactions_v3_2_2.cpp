

// replacement switched on


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

#define COMP_NAMES "names_322.txt" // names of compouds are in this file
#define REP_DG "replace_dG_322.dat"  // dG's of compounds that should be replaced
//#define CONT  "_ext_80_v1"  // suffix added to file names with contributions of groups
#define CONT  "_ext_322"  // suffix added to file names with contributions of groups

// BW #define NUM "v3_ext_rep_80_v2"  // suffix added to all output file names
#define NUM "_4C_v3_2_2_ext_100"  // suffix added to all output file names

void err(char *txt)
{
  cout <<endl<<"err "<<txt<<endl ;
  system("pause") ; exit(0) ;
}

void err(char *txt, string txt2)
{
  cout <<endl<<"err "<<txt<<txt2<<endl ;
  system("pause") ; exit(0) ;
}


void err(char *txt, int c)
{
  cout <<endl<<"err "<<txt<<c<<endl ;
  system("pause") ; exit(0) ;
}


#define SQR(x)  (x)*(x)
long long int xdr48=0x000100010001LL, mndr48=0x0005deece66dLL, doddr48=0xbLL ;
double drand48(void)  // works only on compilers with long long int!
{
  xdr48=mndr48*xdr48+doddr48 ; xdr48&=0xffffffffffffLL ;
  return (xdr48/281474976710656.0) ;
}

void srand48(unsigned int x0)
{
	xdr48=(x0<<16)+0x330e ; 
}

int copy(char *src, char *dest) {
  int i=0 ; do { dest[i]=src[i] ; } while (dest[i++]!=0) ;
  return i-1 ;
}

int copy(const char *src, char *dest) {
  int i=0 ; do { dest[i]=src[i] ; } while (dest[i++]!=0) ;
  return i-1 ;
}

int cmp(char *x, char *y) {
  int i=0 ; 
  do { 
    if (x[i]!=y[i]) break ; 
    i++ ; 
  } while (x[i]!=0 && y[i]!=0) ;
  if (x[i]==0 && y[i]==0) return i ;
  else return 0 ;
}
int cmp(char *x, const char *y) {
  int i=0 ; 
  do { 
    if (x[i]!=y[i]) break ; 
    i++ ; 
  } while (x[i]!=0 && y[i]!=0) ;
  if (x[i]==0 && y[i]==0) return i ;
  else return 0 ;
}

int find_in(char *src, char *txt) 
{
  int i=0,j=0;
  do {
    if (txt[j]==src[i]) {
      j++ ; if (txt[j]==0) return 1 ;
    } else j=0 ;
  } while (src[i++]!=0) ; 
  return 0 ;
}


int _no_groups=0, _e_index=0 ; //, _no_synon=0 ;

struct Atoms {
  int nC,nH,nO,nP,nN ; // nP= no. of phosphorus atoms
  Atoms(int c,int h, int o, int p, int n) { nC=c ; nH=h ; nO=o ; nP=p ; nN=n ; }
  Atoms() { nC=nH=nO=nP=nN=0 ; }
  void zero() { nC=nH=nO=nP=nN=0 ; }
  int mass() { return nC*12+nH+nO*16+31*nP+14*nN ; }

	inline bool operator== (const Atoms& a) const {
		return (nC == a.nC && nH == a.nH && nO == a.nO && nP == a.nP && nN==a.nN);
	}
	inline void operator+= ( const Atoms& a ) {
    nC+=a.nC ; nH+=a.nH ; nO+=a.nO ; nP+=a.nP ; nN+=a.nN ;
  }
	inline void operator-= ( const Atoms& a ) {
    nC-=a.nC ; nH-=a.nH ; nO-=a.nO ; nP-=a.nP ; nN-=a.nN ;
  }
  inline bool operator!= (const int n) const {
		return (nC!=n || nO!=n || nH!=n || nP!=n || nN!=n);     
  }
  inline bool operator!= (const Atoms& a) const {
		return (nC!=a.nC || nO!=a.nO || nH!=a.nH || nP!=a.nP || nN!=a.nN);     
  }
  inline Atoms operator+ (const Atoms& a) const {
    return Atoms(nC+a.nC,nH+a.nH,nO+a.nO,nP+a.nP,nN+a.nN) ;
  }
  inline Atoms operator- (const Atoms& a) const {
    return Atoms(nC-a.nC,nH-a.nH,nO-a.nO,nP-a.nP,nN-a.nN) ;
  }
};  

void cout_formula(Atoms a)
{
  cout <<"C"<<a.nC<<"H"<<a.nH<<"N"<<a.nN<<"O"<<a.nO<<"P"<<a.nP<<" " ;
}

struct Group {
  char name[16] ;   
  vector <string> ligands ; 
  Atoms atoms ;
  int bonds ;
  int mass ;
  int charge ; 
  int i1,i2 ; // indices to the E1/E2 tables, i1 if the group is at the end of the molecule, i2 if the group is in the chain
  
  Group(char *def, char *g1, char *g2=NULL, char *g3=NULL) {
    int i=0,j=0 ;
    ligands.push_back(string(g1)) ;
    if (g2!=NULL) ligands.push_back(string(g2)) ;
    if (g3!=NULL) ligands.push_back(string(g3)) ;
    bonds=0 ; atoms.zero() ; charge=0 ; i1=i2=-1 ;
    do {
      if (def[i]!='-' && def[i]!='=' && def[i]!='#') { 
        name[j++]=def[i] ;
        if (def[i]=='C') if (def[i+1]<'0' || def[i+1]>'9') atoms.nC++ ; else atoms.nC+=def[i+1]-'0' ;
        if (def[i]=='H') if (def[i+1]<'0' || def[i+1]>'9') atoms.nH++ ; else atoms.nH+=def[i+1]-'0' ;
        if (def[i]=='O') if (def[i+1]<'0' || def[i+1]>'9') atoms.nO++ ; else atoms.nO+=def[i+1]-'0' ;
        if (def[i]=='N') if (def[i+1]<'0' || def[i+1]>'9') atoms.nN++ ; else atoms.nN+=def[i+1]-'0' ;
        if (def[i]=='p') { atoms.nP++ ; atoms.nO+=4 ; charge-=2; } // phosphate group -PO_4^{2-}
      } 
      if (def[i]=='-') bonds++ ;
      if (def[i]=='=') bonds+=2 ;
      if (def[i]=='#') bonds+=3 ;
    } while (def[i++]!=0) ; 
    mass=atoms.mass() ; 
    if (bonds==1) i1=_e_index++ ;
    if (bonds==2) { i1=_e_index++ ; i2=_e_index++ ; }
    if (bonds==3) i2=_e_index++ ;
    _no_groups++ ;
  }
  
  int has_ligand(string &txt) {
    for (int i=0;i<ligands.size();i++) if (ligands[i]==txt) return i ;
    return -1 ;
  }
};

Group *gr[]={ new Group("-CH3", "H","H","H"), 
              new Group("-CH2(OH)", "H","H","OH"), // unstable: new Group("-CH(OH)(OH)"), // new Group("-C(OH)(OH)(OH)"),
              new Group("-COOH", "O","OH"), 
              new Group("-CHO", "H","O"),
              new Group("-CH2p", "H","H","p"), 
              new Group("-COp", "O","p"),
              new Group("-CH2(NH2)", "H","H","NH2"), 
              new Group("-CO(NH2)","O","NH2"),
          

              new Group("-CH2-", "H","H"), 
              new Group("-CH(OH)-", "H","OH"),// unstable: new Group("-C(OH)(OH)-"), 
              new Group("-CO-", "O"), 
              new Group("-CHp-", "H","p"),
              new Group("-CH(NH2)-", "H","NH2"),

              new Group("-CH=", "H"), 
              new Group("-C(OH)=", "OH"),
              new Group("-Cp=", "p")
   //           new Group("-C(NH2)=", "NH2")
              
              
//              new Group("-C#") // unlikely

//              new Group("H2CO3")
              
//              new Group("-OP"), new Group("") ; 
                    } ;

// group contributions are loaded into these variables:
double E0,*E1,**E2 ;

void save_groups()
{
  ofstream data("groups_" NUM ".dat") ;
  for (int i=0;i<_no_groups;i++) {
    if (gr[i]->bonds==1) { data <<"-"<<gr[i]->name<<"\t"<<gr[i]->name<<"-"<<endl ; }
    if (gr[i]->bonds==2) { data <<"="<<gr[i]->name<<"\t"<<gr[i]->name<<"="<<endl ; data <<"-"<<gr[i]->name<<"-"<<endl ; }
    if (gr[i]->bonds==3) { data <<"-"<<gr[i]->name<<"="<<"\t"<<"="<<gr[i]->name<<"-"<<endl ; }
//    data<< ;
  }
  data.close() ; 
}         

void load_E_contrib()
{
  double tmp;
  int i,j,no ;
  ifstream data("contrib1" CONT ".dat") ;
  data >> no ; if (no!=_e_index) err("no!=_e_index") ;
  E1=new double[no] ;
  data >> E0 ;
  for (i=0;i<no;i++) data>>E1[i] ;  
  data.close() ; 
  data.clear() ;
  data.open("contrib2" CONT ".dat") ;
  data >> tmp ; E0+=tmp ; //cout <<tmp ; err("x") ;
  for (i=0;i<no;i++) data >> tmp ;
  E2=new double*[no] ;
  for (j=0;j<no;j++) {
    E2[j]=new double[no] ;
    data >> tmp ; for (i=0;i<no;i++) data >> E2[j][i] ;
  }
  data.close() ;
}
                    
Group* find_group(char *txt)
{
  for (int i=0;i<_no_groups;i++) if (cmp(gr[i]->name,txt)) return gr[i] ;
  cout <<endl<<"find_group: cannot find "<<txt<<endl ; system("pause") ; exit(0) ;
}                    

/*Group* remove_hydrogen(Group *x)
{
  char *txt = x->string ;
  if (cmp(txt,"CH3")) return find_group("CH2") ;
  if (cmp(txt,"CH2p")) return find_group("CHp") ;
  if (cmp(txt,"CH(OH)(OH)")) return find_group("C(OH)(OH)") ;
  if (cmp(txt,"CH2(OH)")) return find_group("CH(OH)") ;
  if (cmp(txt,"CH2")) return find_group("CH") ;
  if (cmp(txt,"CHp")) return find_group("Cp") ;
  if (cmp(txt,"CH(OH)")) return find_group("C(OH)") ;
  
  cout <<endl<<"remove_hydrogen: cannot find "<<txt<<endl ; system("pause") ; exit(0) ;
}*/

Group* add_hydrogen_to_C(Group *x, int show_error=1)
{
  // this should be modified to scan through ligands and find automatically the correct one
//  Group *tmp=NULL ;
/*  int i;
  for (i=0;i<_no_groups;i++) if (x->atoms+Atoms(0,1,0,0,0)==gr[i]->atoms) return gr[i] ;
  if (show_error) { cout <<endl<<"add_hydrogen: cannot find "<<x->name<<endl ; system("pause") ; exit(0) ; }
  else return 0 ;
  */

// for tests:
//  for (int i=0;i<_no_groups;i++) if (x->atoms+Atoms(0,1,0,0,0)==gr[i]->atoms) { tmp= gr[i] ; break ; }
//  if (show_error) { cout <<endl<<"add_hydrogen: cannot find "<<x->name<<endl ; system("pause") ; exit(0) ; }
  
  int i,j,ok;
  for (i=0;i<_no_groups;i++) {
    if (x->ligands.size()==gr[i]->ligands.size()-1) {
      ok=1 ;
      for (j=0;j<x->ligands.size();j++) if (x->ligands[j]!=gr[i]->ligands[j+1]) ok=0 ; 
      if (ok && gr[i]->ligands[0]=="H") return gr[i] ;
    }
  }
  if (show_error) { cout <<endl<<"add_hydrogen: cannot find "<<x->name<<endl ; system("pause") ; exit(0) ; }
  else return 0 ; // { if (tmp!=NULL) { cout <<x->name<<"-->"<<tmp->name ; err("y") ; } return 0 ; }
}

/*Group* add_hydrogen(Group *x, int show_error=1)
{
  char *txt = x->string ;
  if (cmp(txt,"CH2")) return find_group("CH3") ;
  if (cmp(txt,"C(OH)(OH)")) return find_group("CH(OH)(OH)") ;
  if (cmp(txt,"CH(OH)")) return find_group("CH2(OH)") ;
  if (cmp(txt,"CH")) return find_group("CH2") ;
  if (cmp(txt,"C(OH)")) return find_group("CH(OH)") ;
  if (cmp(txt,"CO")) return find_group("CHO") ;
  
  if (show_error) { cout <<endl<<"add_hydrogen: cannot find "<<txt<<endl ; system("pause") ; exit(0) ; }
  else return 0 ;
}*/

                    
/*Group* remove_hydroxyl(Group *x)
{
//  int i;
//  for (i=0;i<_no_groups;i++) if (x==gr[i]) break ;
//  if (i>=_no_groups) err("cannot find the source group") ; 
  char *txt = x->string ;
  if (cmp(txt,"CH2(OH)")) return find_group("CH2") ;
  if (cmp(txt,"CH(OH)(OH)")) return find_group("CH(OH)") ;
  if (cmp(txt,"CH(OH)")) return find_group("CH") ;
  if (cmp(txt,"C(OH)(OH)")) return find_group("C(OH)") ;
  
  cout <<endl<<"remove_hydroxyl: cannot find "<<txt<<endl ; system("pause") ; exit(0) ;
}*/

/*class Synonyms {
  public:
    Group *src1, *src2, *dest1, *dest2 ;
    Synonyms::Synonyms(char *g1, char *g2, char *g3, char *g4=NULL) {
      int i ;
      src1=src2=dest1=dest2=NULL ;
      for (i=0;i<_no_groups;i++) {
        if (cmp(gr[i]->string,g1)) src1=gr[i] ;
        if (cmp(gr[i]->string,g2)) src2=gr[i] ;
        if (cmp(gr[i]->string,g3)) dest1=gr[i] ;
        if (g4!=NULL && cmp(gr[i]->string,g4)) dest2=gr[i] ;
      }
      if (src1==NULL || src2==NULL || (dest1==NULL && dest2==NULL)) err("synonyms 1") ;
      _no_synon++ ;
    }
};


Synonyms *synonyms[]={new Synonyms("C(COOH)","OH","C(OH)(COOH)"),new Synonyms("C(OH)","COOH","C(OH)(COOH)"),
                      new Synonyms("CH(COOH)","OH", "CH(OH)","COOH")} ;
*/
class Compound {
  public:
    char string[256] ; // string representation (molecular formula) of the compound
    char name[256] ; // name of the molecule
    Group *groups[32] ; // group representation of the compound
    int free_bonds[32] ; // no. of free bonds for each group
    int no_groups ;
    //int left_bonds, right_bonds ;
    Atoms atoms ;
    int mass ;
    int charge ; // charge 
    double dG ; // gibbs energy of formation \Delta G^0_f
        
    Compound() {
      no_groups=mass=charge=0 ; atoms.zero() ; dG=0 ; name[0]=0 ;
    }
    void clear()
    {
      string[0]=0 ; no_groups=mass=charge=0 ; atoms.zero() ; dG=0 ; name[0]=0 ;
    }
    Compound(char str[], double deltag, int ch, int c,int h,int o,int p, int n) {
      charge=ch ; atoms.nC=c ; atoms.nH=h ; atoms.nO=o ; atoms.nP=p ; atoms.nN=n ; mass=atoms.mass() ;
      no_groups=0 ; dG=deltag ; name[0]=0 ;
      copy(str,string) ;
    }
    Compound* clone() {
      Compound *cp=new Compound ;
      copy(string,cp->string) ;
      cp->no_groups=no_groups ;
      for (int i=0;i<no_groups;i++) { cp->groups[i]=groups[i] ; cp->free_bonds[i]=free_bonds[i] ; }
//      cp->left_bonds=left_bonds ; cp->right_bonds=right_bonds ;
      cp->atoms=atoms ;
      cp->mass=mass ; cp->dG=dG ; cp->charge=charge ; copy(name,cp->name) ;
      return cp ;     
    }    
    void add_name() {
      char str[256],txt[256],txt2[256] ;
      FILE *names=fopen(COMP_NAMES,"r") ;
      do {
        fscanf(names,"%s %s %s",str,txt,txt2) ; // txt2= second name corresponding to a function in the mathematica notebook
        if (cmp(str,string)) { copy (txt,name) ; break ; }
      } while (!feof(names)) ;
      fclose(names);
    }
    void if_possible_replace_dG_by_Alberty() {
      char str[256] ; float deltag ;
      FILE *names=fopen(REP_DG,"r") ; if (names==NULL) err("file " REP_DG " not found! ") ;
      do {
        fscanf(names,"%s %f",str,&deltag) ; 
        if (cmp(str,string)) { dG=deltag ; break ; }
      } while (!feof(names)) ;
      fclose(names);
      
    }
    void update() {
      mass=charge=0 ; atoms.zero() ; 
      for (int i=0;i<no_groups;i++) {
        atoms+=groups[i]->atoms ; charge+=groups[i]->charge ;
      }
      mass=atoms.mass() ;
    }
    int dangling_bonds() {
      int sum=0 ; for (int i=0;i<no_groups;i++) sum+=free_bonds[i] ;
      return sum ;
    }
/*    Compound::Compound(char *def) {
      copy(def,string) ; this->find_mol_formula() ;
    }*/
    int add_right(Group *g) {
//      cout <<g->string<<"\t" ;
      if (no_groups==0) { 
        free_bonds[no_groups]=g->bonds ;
        groups[no_groups++]=g ; 
//        if (g->bonds>=1) right_bonds=1 ; 
//        if (g->bonds==2) left_bonds=1 ;
        atoms+=g->atoms ;
        mass=atoms.mass() ; charge+=g->charge ;
        return 0 ; 
      }
      else if (free_bonds[no_groups-1]>0 && g->bonds>=free_bonds[no_groups-1]) {
//        int connected = MIN(free_bonds[no_groups-1], g->bonds) ;
        free_bonds[no_groups]=g->bonds - free_bonds[no_groups-1] ;
        free_bonds[no_groups-1]=0 ;
        groups[no_groups++]=g ;
        //right_bonds=g->bonds-1 ;
        atoms+=g->atoms ; 
        mass=atoms.mass() ;  charge+=g->charge ;
      }
      return -1 ;
    }
    void remove_right() {
      if (no_groups==0) return ;
      Group *g=groups[no_groups-1] ;
      atoms-=g->atoms ; 
      mass=atoms.mass() ; charge-=g->charge ;
      if (no_groups>1) free_bonds[no_groups-2]=g->bonds-free_bonds[no_groups-1] ;
      else free_bonds[0]=0 ;
      no_groups-- ;
    }
    double free_energy() {
//      int prev=0 ;
      dG=E0 ;//53.88+groups[0]->dG1 ;
      dG+=E1[groups[0]->i1]  ;
      for (int i=1;i<no_groups-1;i++) dG+=E1[groups[i]->i2] ;
      dG+=E1[groups[no_groups-1]->i1] ;
      
      if (no_groups==2) dG+=E2[groups[0]->i1][groups[1]->i1] ;
      else if (no_groups==3) {
        if (groups[0]->i1==-1) err("0, i1") ;
        if (groups[1]->i2==-1) { cout <<groups[1]->name ; err("1, i2") ; }
        if (groups[2]->i1==-1) err("2, i1") ;
        dG+=E2[groups[0]->i1][groups[1]->i2] + E2[groups[1]->i2][groups[2]->i1] ; 
      } else { 
        if (groups[0]->i1==-1) err("0, i1") ;
        for (int i=1;i<no_groups-1;i++) if (groups[1]->i2==-1) { cout <<groups[1]->name ; err("i, i2, i=",i) ; }
        if (groups[no_groups-1]->i1==-1) err("no_groups-1, i1") ;
        dG+=E2[groups[0]->i1][groups[1]->i2] + E2[groups[no_groups-2]->i2][groups[no_groups-1]->i1] ;
        for (int i=1;i<no_groups-2;i++) dG+=E2[groups[i]->i2][groups[i+1]->i2] ;
        //err("more than 3 groups, dG cannot be calculated") ;
      }
//      cout<<endl ;       
    }
    int charge_in_solution() {
      int c=charge ;
      for (int i=0;i<no_groups;i++) if (groups[i]==find_group("COOH")) c-=1 ;
      return c ; 
    }
    void optimize(int err_if_dangling=1) {
//      this->update() ;
      if (no_groups<2) return ;

      for (int i=0;i<no_groups;i++) if (free_bonds[i]>0) if (err_if_dangling) err("cannot optimize with tangling bonds (i=",i) ; else return ;

/*      for (int j=0;j<_no_synon;j++) {
        for (int i=0;i<no_groups-1;i++) {
          if ((groups[i]==synonyms[j]->src1 && groups[i+1]==synonyms[j]->src2) || 
            (groups[i]==synonyms[j]->src2 && groups[i+1]==synonyms[j]->src1)) {
              if (synonyms[j]->dest2==NULL) { // replace two groups by one group  
                groups[i]=synonyms[j]->dest1 ; 
                for (int k=i+1;k<no_groups-1;k++) groups[k]=groups[k+1] ;
                no_groups-- ;  //err("x") ;
              } else { // replace two groups by two other groups
                if (groups[i]==synonyms[j]->src1) { groups[i]=synonyms[j]->dest1 ; groups[i+1]=synonyms[j]->dest2 ; }
                else { groups[i+1]=synonyms[j]->dest1 ; groups[i]=synonyms[j]->dest2 ; }
              }
            }
        }
      }     */ 
      
      for (int j=0;j<no_groups/2;j++) {
        if (groups[j]->mass > groups[no_groups-1-j]->mass) { // revert the sequence
          for (int i=0;i<no_groups/2;i++) { 
            Group *tmp=groups[i] ; int tmpf=free_bonds[i] ;
            groups[i]=groups[no_groups-1-i] ; groups[no_groups-1-i]=tmp ;
            free_bonds[i]=free_bonds[no_groups-1-i] ; free_bonds[no_groups-1-i]=tmpf ;
          }
          break ;
        } else if (groups[j]->mass < groups[no_groups-1-j]->mass) break ;
      }  
    }
    void to_string() {
      if (no_groups==0) { string[0]=0 ; return ; }
      for (int i=0;i<no_groups-1;i++) if (free_bonds[i]>0) err("cannot generate string with dangling bonds inside the molecule") ;
      int j=0 ;
//      cout <<"x"<<no_groups<<" " ;
//      if (left_bonds) string[j++]='-' ; 
      int prev=groups[0]->bonds ;
      for (int i=0;i<no_groups;i++) {
        j+=copy(groups[i]->name,string+j) ; 
        if (i<no_groups-1 && prev==1) string[j++]='-' ; 
        if (i<no_groups-1 && prev==2) string[j++]='=' ; 
        if (i<no_groups-1 && prev==3) string[j++]='#' ; 
        if (i<no_groups-1 && prev==0) err("? this should never happen") ; //string[j++]='?' ;
        if (i<no_groups-1) prev=groups[i+1]->bonds-prev ;
      } 
//      string[j++]=0 ; 
      if (free_bonds[no_groups-1]>0) {
        if (free_bonds[no_groups-1]==1) string[j++]='-' ;
        if (free_bonds[no_groups-1]==2) string[j++]='=' ;
        if (free_bonds[no_groups-1]==3) string[j++]='#' ;
        string[j++]=0 ;
      }
      //if (right_bonds) { string[j++]='-' ; string[j++]=0 ; } 
    }
    void print() {
      this->to_string() ;
      cout <<string<<endl ;
//      for (int i=0;i<no_groups;i++) cout <<groups[i]->string<<" "; cout <<endl ; 
    }
    
    // this function should not be used except for constructing reactions because it does not update other elements of Compound
    void add_hydrogen_to_last() {
      Group *g=add_hydrogen_to_C(groups[no_groups-1],0) ; 
      if (g!=NULL) { remove_right() ; add_right(g) ; } //err("x") ; }
    }
};

// char name, double deltaG,  int ch, int c,int h,int o,int p, int n
// values of deltaG copied from a Mathematica notebook, valid for ph=7 and is=0.2
Compound *_water=new Compound("H2O",-155.758, 0, 0,2,1,0,0) ;
Compound *_co2=new Compound("CO2(aq)",-547.111, 0, +1,2,3,0,0) ;  // this is indeed H2CO3
Compound *_nad=new Compound("NAD_ox",1057.86, +1, 0,0,0,0,0) ; // this is NAD+, oxidised version of NAD
Compound *_nadh=new Compound("NAD_red",1118.95, +1, 0,2,0,0,0) ; // this NADH + H+, but this is at equilibrium with NADH2+ 
                                                                // so protons are automatically included 
Compound *_amp=new Compound("AMP",-555.258, 0, 0,2,4,1,0) ; // this is AMP + 2H+
Compound *_adp=new Compound("ADP",-1424.91, -2, 0,1,7,2,0) ; // this is ADP + H+
Compound *_atp=new Compound("ATP",-2292.39, -4, 0,0,10,3,0) ;
Compound *_pi=new Compound("Pi", -1059.4, -2, 0,1,4,1,0) ;  // this is PO_4^{3-} + H+
Compound *_ppi=new Compound("PPi", -1940.27, -4, 0,0,7,2,0) ;  // this is P_2O_7^{4-}    

Compound *_nh3=new Compound("NH3(aq)",82.78, 0, 0,3,0,0,1) ; // this is in fact NH3(aq) but represented as NH3
Compound *_urea=new Compound("CO(NH2)2",-39.93, 0, 1,4,1,0,2) ; // urea

Compound *_nh2donor=new Compound("NH2donor",-372.5, -2, 5,9,4,0,1) ; // for example glutamate
Compound *_nh2acceptor=new Compound("NH2acceptor",-633.585, -2, 5,6,5,0,0) ; // for example 2-Oxoglutarate (ketoglutarate)
 

Compound *cc_basic[13]={_water,_co2,_nad,_nadh,_amp,_adp,_atp,_pi,_ppi,_nh3,_urea,_nh2donor,_nh2acceptor} ;
int _no_cc_basic=13 ;

int equalQ(Compound *x, Compound *y) 
{
  if (x->no_groups!=y->no_groups) return 0 ;
  for (int i=0;i<x->no_groups;i++) if (x->groups[i]!=y->groups[i]) return 0 ;
  return 1 ; 
}

int isomersQ(Compound *x, Compound *y) 
// returns 1 if x and y are isomers, i.e., they have the same molecular formula
// but different structural formulas,
// otherwise it returns 0
{
  if (x->atoms==y->atoms) return 1 ;
  else return 0 ;
}


class Reaction {
  public:
    Compound *substr[4], *prod[4] ;
    int no_substr, no_prod ;
    char string[256], name[64] ;
    double dG ;
    Reaction(char *txt) {
      no_substr=no_prod=0 ; dG=0 ;
      copy(txt,name) ;
    }
    Reaction(const char *txt) {
      no_substr=no_prod=0 ; dG=0 ;
      copy(txt,name) ;
    }
    void substrate(Compound *c) { substr[no_substr++]=c ; }
    void product(Compound *c) { prod[no_prod++]=c ; }
    void to_string() {
      if (no_substr==1 && no_prod==1) sprintf(string,"%s ---> %s",substr[0]->string,prod[0]->string) ; 
      if (no_substr==1 && no_prod==2) sprintf(string,"%s ---> %s + %s",substr[0]->string,prod[0]->string,prod[1]->string) ; 
      if (no_substr==2 && no_prod==1) sprintf(string,"%s + %s ---> %s",substr[0]->string,substr[1]->string,prod[0]->string) ; 
      if (no_substr==2 && no_prod==2) sprintf(string,"%s + %s ---> %s + %s",substr[0]->string,substr[1]->string,prod[0]->string,prod[1]->string) ; 
      if (no_substr==3 && no_prod==2) sprintf(string,"%s + %s + %s ---> %s + %s",substr[0]->string,substr[1]->string,substr[2]->string,prod[0]->string,prod[1]->string) ; 
      if (no_substr==2 && no_prod==3) sprintf(string,"%s + %s ---> %s + %s + %s",substr[0]->string,substr[1]->string,prod[0]->string,prod[1]->string,prod[2]->string) ; 
      if (no_substr==3 && no_prod==3) sprintf(string,"%s + %s + %s ---> %s + %s + %s",substr[0]->string,substr[1]->string,substr[2]->string,prod[0]->string,prod[1]->string,prod[2]->string) ; 
      if (no_substr==4 && no_prod==3) sprintf(string,"%s + %s + %s + %s ---> %s + %s + %s",substr[0]->string,substr[1]->string,substr[2]->string,substr[3]->string,prod[0]->string,prod[1]->string,prod[2]->string) ; 
      if (no_substr==3 && no_prod==4) sprintf(string,"%s + %s + %s ---> %s + %s + %s + %s",substr[0]->string,substr[1]->string,substr[2]->string,prod[0]->string,prod[1]->string,prod[2]->string,prod[3]->string) ; 
      if (no_substr==4 && no_prod==4) sprintf(string,"%s + %s + %s + %s ---> %s + %s + %s + %s",substr[0]->string,substr[1]->string,substr[2]->string,substr[3]->string,prod[0]->string,prod[1]->string,prod[2]->string,prod[3]->string) ; 
      if (no_substr>4 || no_prod>4) err("too many substrates/products") ;
    }
    void update() {
      dG=0 ; for (int i=0;i<no_prod;i++) dG+=prod[i]->dG ; 
      for (int i=0;i<no_substr;i++) dG-=substr[i]->dG ;
      this->to_string() ; 
  
//      return ;    
      // check if all atoms and charges are balanced
      int charge=0, mass=0 ;
      Atoms atoms ;
      for (int i=0;i<no_prod;i++) { atoms+=prod[i]->atoms ; charge+=prod[i]->charge ; mass+=prod[i]->mass ; }
      for (int i=0;i<no_substr;i++) { atoms-=substr[i]->atoms ; charge-=substr[i]->charge ; mass-=substr[i]->mass ; }
      if (charge!=0 || atoms!=0 || mass!=0) {
        cout <<"\nerror in reaction:\n"<<string<<"\n dmass="<<mass<<" dcharge="<<charge<<" dnC="<<atoms.nC<<" dnH="<<atoms.nH<<" dnO="<<atoms.nO<<" dnP="<<atoms.nP<<endl ;
        for (int i=0;i<no_substr;i++) cout <<substr[i]->charge<<" " ; cout <<"--->" ;
        for (int i=0;i<no_prod;i++) cout <<prod[i]->charge<<" " ; cout <<endl ;
        system("pause");
        exit(1) ;
      }
    }
};

template<typename T>
void reverse(vector<T> *a)
{
  T tmp ; 
  int n=a->size() ;
  for (int i=0;i<n/2;i++) {
//    cout <<"*"<<n-1-i ;
    tmp=(*a)[i] ; (*a)[i]=(*a)[n-1-i] ; (*a)[n-1-i]=tmp ;
  }
//  cout <<"--"<<(*a)[1]<<"--" ;
}

int parse(string s, vector <string> *v, char *sep) 
{
  int i,j=0,k=0;
  for (i=0;i<s.length();i++) {
    j=0 ; while (sep[j]!=0) {
      if (s[i]==sep[j]) { v->push_back(s.substr(k,i-k)) ; k=i+1 ; } 
      j++ ;
    }
  }
  if (k>0) v->push_back(s.substr(k,i-k)) ;
  return v->size() ;
}



struct Reaction_class {
  string human_readable ; 
  vector <Compound*> ext_p, ext_s ; // external metabolites
  vector <int> gr_p, gr_s ; // groups of both compounds, -1,-2,... if R1,R2,...
  vector <char> b_p, b_s ; // bonds
  vector <string> string_p, string_s ; // strings denoting groups. Important for R1,R2
  vector <int> grat_p, grat_s ; // group attributes: 1 = at the end
  Reaction_class() {
    human_readable="" ; 
  }
  Reaction_class(Reaction_class *src) {
    human_readable=src->human_readable ;
    ext_p=src->ext_p ; ext_s=src->ext_s ;
    gr_p=src->gr_p ; gr_s=src->gr_s ;
    b_p=src->b_p ; b_s=src->b_s ;
    string_p=src->string_p ; string_s=src->string_s ;
    grat_p=src->grat_p ; grat_s=src->grat_s ;
  }
  Reaction_class(string s) {
    human_readable=s ;
//    cout <<"new reaction: "<<s<<endl<<endl ;
    ext_p.clear() ; ext_s.clear() ;
    vector <string> parsed ;
    parse(s,&parsed," ") ;
    //char *parsed=strtok(s," ") ;
    int i,j,k,side=0;
    for (i=0;i<parsed.size();i++) {
//      cout <<parsed[i]<<"_   " ; 
      if (parsed[i]==">") { side=1 ; /*cout <<endl ; */continue ; }
      for (j=0;j<_no_cc_basic;j++) if (cmp(cc_basic[j]->string,parsed[i].data())) break ;
      if (j<_no_cc_basic) { /*cout <<"ext"<<j ; */if (side==0) ext_s.push_back(cc_basic[j]) ; else ext_p.push_back(cc_basic[j]) ; }
      else {
        vector <string> grs ;
        parse(parsed[i],&grs,"-=#") ; 
        if (grs.size()==0 && parsed[i].length()>=2) grs.push_back(parsed[i]) ;
        if (grs.size()>1) { // if at least two groups then determine the number of bonds
          for (j=0;j<parsed[i].length();j++) {
            if (parsed[i][j]=='-' || parsed[i][j]=='=' || parsed[i][j]=='#') { 
              if (side==0) b_s.push_back(parsed[i][j]) ; else b_p.push_back(parsed[i][j]) ;
            }
          }
        }
        for (j=0;j<grs.size();j++) {
//          cout <<grs[j]<<" " ;
          for (k=0;k<_no_groups;k++) if (cmp(gr[k]->name,grs[j].data())) break ;
          if (k<_no_groups) { // fixed groups
            int attr=0 ; 
            if (j==0 || j==grs.size()-1) attr=1 ; // group at the end
//            if (attr && gr[k]->bonds!=1) err("eeee") ;
            if (!attr && gr[k]->bonds<=1) { human_readable="" ; return ; } // exit if end groups in the middle - used when constructing new reactions from wildcard patterns
            if (side==0) { 
              gr_s.push_back(k) ; grat_s.push_back(attr) ; string_s.push_back(gr[k]->name) ;
            } else { 
              gr_p.push_back(k) ; grat_p.push_back(attr) ; string_p.push_back(gr[k]->name) ;
            }
          } else { // R1, R2, ...
            if (grs[j][0]=='R') {
              if (grs[j][1]<49 || grs[j][1]>57) err("R is missing an index") ;
              if (side==0) {
                gr_s.push_back(48-grs[j][1]) ; grat_s.push_back(1) ; string_s.push_back(grs[j]) ; // gr_s[] = -1,-2,...
              } else {
                gr_p.push_back(48-grs[j][1]) ;  grat_p.push_back(1) ; string_p.push_back(grs[j]) ; // gr_p[] = -1,-2,...
              }
            } else err("unidentified group:", grs[j]) ;
          }
        }
      }      
//      cout<<endl ;
    }
//    for (i=0;i<gr_s.size();i++) if (gr_s[i]>=0) cout <<gr[gr_s[i]]->name<<" " ; else cout <<gr_s[i]<<" " ; cout <<endl ;
//    for (i=0;i<gr_p.size();i++) if (gr_p[i]>=0) cout <<gr[gr_p[i]]->name<<" " ; else cout <<gr_p[i]<<" " ; cout <<endl ;
//    err("x",gr_s.size()) ;
    //parse_compound(s) 
  }
  int check_if_reaction_possible(Compound *x, Compound *y) ;
  void update_human_readable() ;
  bool invalid() ;
};

void Reaction_class::update_human_readable() 
{
  human_readable="" ;
  for (int i=0;i<gr_s.size();i++) {
    if (gr_s[i]>=0) human_readable+=gr[gr_s[i]]->name ; else human_readable+=string_s[i] ; 
    if (i<b_s.size()) human_readable+=b_s[i] ;
  }
  for (int i=0;i<ext_s.size();i++) human_readable+=" + "+string(ext_s[i]->string) ;
  human_readable+=" > " ;
  for (int i=0;i<gr_p.size();i++) {
    if (gr_p[i]>=0) human_readable+=gr[gr_p[i]]->name ; else human_readable+=string_p[i] ; 
    if (i<b_p.size()) human_readable+=b_p[i] ;
  }    
  for (int i=0;i<ext_p.size();i++) human_readable+=" + "+string(ext_p[i]->string) ;
}

bool Reaction_class::invalid()
{
  cout <<"rr "<<human_readable<<endl ;
  if (human_readable=="") return true ;
  if (gr_s.size()!=b_s.size()+1) err("invalid(): dangling bonds") ;  
  if (gr_p.size()!=b_p.size()+1) err("invalid(): dangling bonds") ;  
  for (int i=0;i<gr_s.size();i++) {
    if (gr_s[i]>=0) {
      int b=0 ;
      if (i<gr_s.size()-1) b+=(b_s[i]=='-'?1:0)+(b_s[i]=='='?2:0)+(b_s[i]=='#'?3:0) ;
      if (i>0) b+=(b_s[i-1]=='-'?1:0)+(b_s[i-1]=='='?2:0)+(b_s[i-1]=='#'?3:0) ;
      if (b!=gr[gr_s[i]]->bonds) return true ;
    }
  }
  for (int i=0;i<gr_p.size();i++) {
    if (gr_p[i]>=0) {
      int b=0 ;
      if (i<gr_p.size()-1) b+=(b_p[i]=='-'?1:0)+(b_p[i]=='='?2:0)+(b_p[i]=='#'?3:0) ;
      if (i>0) b+=(b_p[i-1]=='-'?1:0)+(b_p[i-1]=='='?2:0)+(b_p[i-1]=='#'?3:0) ;
      if (b!=gr[gr_p[i]]->bonds) return true ;
    }
  }
  return false ;
}


struct Reaction_classes {
  vector <Reaction_class*> rc ;
  string name ;
  Reaction_classes() { name="" ; rc.clear() ; }
  Reaction_classes(string n) { name=n ; rc.clear() ; }
  Reaction_classes(string n, char *s1) { name=n ; rc.push_back(new Reaction_class(s1)) ; }
  Reaction_classes(string n, char *s1, char *s2) { name=n ; rc.push_back(new Reaction_class(s1)) ; rc.push_back(new Reaction_class(s2)) ; }
  Reaction_classes(string n, char *s1, char *s2, char *s3) { name=n ; rc.push_back(new Reaction_class(s1)) ; rc.push_back(new Reaction_class(s2)) ; rc.push_back(new Reaction_class(s3)) ; }
  Reaction_classes(string n, char *s1, char *s2, char *s3, char *s4) { name=n ; rc.push_back(new Reaction_class(s1)) ; rc.push_back(new Reaction_class(s2)) ; 
    rc.push_back(new Reaction_class(s3)) ; rc.push_back(new Reaction_class(s4)) ; }
  Reaction_classes(string n, char *s1, char *s2, char *s3, char *s4, char *s5) { name=n ; rc.push_back(new Reaction_class(s1)) ; rc.push_back(new Reaction_class(s2)) ; 
    rc.push_back(new Reaction_class(s3)) ; rc.push_back(new Reaction_class(s4)) ; rc.push_back(new Reaction_class(s5)) ; }
  Reaction_classes(string n, char *s1, char *s2, char *s3, char *s4, char *s5, char *s6) { name=n ; rc.push_back(new Reaction_class(s1)) ; rc.push_back(new Reaction_class(s2)) ; 
    rc.push_back(new Reaction_class(s3)) ; rc.push_back(new Reaction_class(s4)) ; rc.push_back(new Reaction_class(s5)) ; rc.push_back(new Reaction_class(s6)) ; }
  void add(const char *s) { rc.push_back(new Reaction_class(s)) ; }  
  void remove_repeated_reactions() {
    for (int i=0;i<rc.size();i++)
      for (int j=i+1;j<rc.size();j++) 
        if (rc[i]->human_readable==rc[j]->human_readable) { rc.erase(rc.begin()+j) ; j-- ; }
  }
  void remove_invalid_reactions() {
    for (int i=0;i<rc.size();i++) if (rc[i]->invalid()) { rc.erase(rc.begin()+i) ; i-- ; }
  }
  void add_more_by_replacing_R_by_H() {
    for (int i=0;i<rc.size();i++) {
      Reaction_class *r=rc[i] ;
//      if (r==NULL) err("null!") ;
      if (r->gr_s.size()>1 && r->b_p.size()>1 && r->b_s.size()>0 && r->b_p.size()>0 
          && r->gr_s[r->gr_s.size()-1]==-2 && r->b_s[r->b_s.size()-1]=='-' && r->b_p[r->b_p.size()-1]=='-'
          && r->gr_p[r->gr_p.size()-1]==-2) { // substr and prod both contain -R2 at the end {
//        r->update_human_readable() ;
        Reaction_class *nr=new Reaction_class(r) ;
        cout <<"old:"<<nr->human_readable<<endl ;
        nr->gr_s.pop_back() ; nr->b_s.pop_back() ; nr->string_s.pop_back() ; nr->grat_s.pop_back() ;
        int n=nr->gr_s.size()-1,g ;
        Group *gs=add_hydrogen_to_C(gr[nr->gr_s[n]],0) ;
        if (gs!=NULL) {
          for (g=0;g<_no_groups;g++) if (gr[g]==gs) break ;
          if (g==_no_groups) err("!!!") ;
          nr->gr_s[n]=g ;
          nr->string_s[n]=gr[g]->name ;
          nr->grat_s[n]=1 ;

          nr->gr_p.pop_back() ; nr->b_p.pop_back() ; nr->string_p.pop_back() ; nr->grat_p.pop_back() ;
          int n=nr->gr_p.size()-1,g ;
          Group *gp=add_hydrogen_to_C(gr[nr->gr_p[n]],0) ;
          if (gp!=NULL) {
            for (g=0;g<_no_groups;g++) if (gr[g]==gp) break ;
            if (g==_no_groups) err("!!!") ;
            nr->gr_p[n]=g ;
            nr->string_p[n]=gr[g]->name ;
            nr->grat_p[n]=1 ;        
            nr->update_human_readable() ;        
            rc.push_back(nr) ; 
            cout <<"new:"<<nr->human_readable<<endl ;
          }
        }
      }
      if (r->gr_s.size()>1 && r->b_s.size()>0 && r->b_p.size()>0 && r->gr_p.size()>1
          && r->gr_s[0]==-1 && r->b_s[0]=='-' && r->b_p[0]=='-' && r->gr_p[0]==-1) { // substr and prod both contain R1- at the end 
//        r->update_human_readable() ;

        Reaction_class *nr=new Reaction_class(r) ;
        cout <<"old2:"<<nr->human_readable<<endl ;
        nr->gr_s.erase(nr->gr_s.begin()) ; nr->b_s.erase(nr->b_s.begin()) ; nr->string_s.erase(nr->string_s.begin()) ; nr->grat_s.erase(nr->grat_s.begin()) ;
        int g ;
        Group *gs=add_hydrogen_to_C(gr[nr->gr_s[0]],0) ;
        if (gs!=NULL) {
          for (g=0;g<_no_groups;g++) if (gr[g]==gs) break ;
          if (g==_no_groups) err("!!!") ;
          nr->gr_s[0]=g ;
          nr->string_s[0]=gr[g]->name ;
          nr->grat_s[0]=1 ;

          nr->gr_p.erase(nr->gr_p.begin()) ; nr->b_p.erase(nr->b_p.begin()) ; nr->string_p.erase(nr->string_p.begin()) ; nr->grat_p.erase(nr->grat_p.begin()) ;
          int g ;
          Group *gp=add_hydrogen_to_C(gr[nr->gr_p[0]],0) ;
          if (gp!=NULL) {
            for (g=0;g<_no_groups;g++) if (gr[g]==gp) break ;
            if (g==_no_groups) err("!!!") ;
            nr->gr_p[0]=g ;
            nr->string_p[0]=gr[g]->name ;
            nr->grat_p[0]=1 ;        
            // now we need to reverse if it contains -R2
            if (nr->gr_s.size()>1 && nr->gr_s[nr->gr_s.size()-1]==-2) {
              reverse(&nr->gr_s) ; reverse(&nr->b_s) ;reverse(&nr->string_s) ; reverse(&nr->grat_s) ;
              reverse(&nr->gr_p) ; reverse(&nr->b_p) ;reverse(&nr->string_p) ; reverse(&nr->grat_p) ;
              nr->gr_s[0]=-1 ; nr->gr_p[0]=-1 ;
              nr->string_s[0][1]='1' ; nr->string_p[0][1]='1' ;
            }
            nr->update_human_readable() ;        
            rc.push_back(nr) ; 
            cout <<"new2:"<<nr->human_readable<<endl ;
          }
        }
      }
    }
  }
};


int Reaction_class::check_if_reaction_possible(Compound *x, Compound *y) 
{
  int ok=0, diff=1, v=0 ;

//  if (human_readable=="R1-CH2(NH2) + H2O + H2O > R1-CH2(OH) + NH3(aq)") v=1 ;

/*  vector <string> ss ;
  ss.push_back("11") ; ss.push_back("22") ;
  for (i=0;i<ss.size();i++) cout <<ss[i] ; cout <<endl ;
  reverse(&ss) ;  
  for (i=0;i<ss.size();i++) cout <<ss[i] ; cout <<endl ;
*/
//  cout <<"x" ;

  // first check if the number of atoms agrees, if not then exit
  Atoms a=y->atoms-x->atoms ; 
//  for (i=0;i<rc->gr_p.size();i++) { nC+=gr[rc->gr_p[i]]->nC ; nH+=gr[rc->gr_p[i]]->nH ; nO+=gr[rc->gr_p[i]]->nO ; nP+=gr[rc->gr_p[i]]->nP ; }
//  for (i=0;i<rc->gr_s.size();i++) { nC-=gr[rc->gr_s[i]]->nC ; nH-=gr[rc->gr_s[i]]->nH ; nO-=gr[rc->gr_s[i]]->nO ; nP-=gr[rc->gr_s[i]]->nP ; }
  for (int i=0;i<ext_p.size();i++) a+=ext_p[i]->atoms ;
  for (int i=0;i<ext_s.size();i++) a-=ext_s[i]->atoms ; 
/*  if (a!=0 && v && x->atoms.nN>0) {
    cout_formula(x->atoms) ; for (int i=0;i<ext_s.size();i++) cout_formula(ext_s[i]->atoms) ;
    cout <<endl ;
    cout_formula(y->atoms) ; for (int i=0;i<ext_p.size();i++) cout_formula(ext_p[i]->atoms) ;
    err("x") ;
  }*/
  if (a!=0) return 0 ;

  Compound Rs[2],Rp[2] ; // Rs

//  cout <<"type:"<<human_readable<<"\n" ;
  if (v) cout <<endl<<"reaction: "<<x->string<<" "<<y->string<<" : " ;


// search for pattern in the substrate  
  for (int s1=0;s1<2;s1++) {
    int n=gr_s.size() ;
    for (int is=0;is<=x->no_groups-n;is++) {
      diff=0 ; //cout <<"£" ;
      for (int j=0;j<n;j++) {
        if (gr_s[j]>=0 && (gr[gr_s[j]]!=x->groups[is+j] || (grat_s[j]==1 && is+j>0 && is+j<x->no_groups-1)/* not at the end*/)) { diff=1 ; break ; }
      }
      if (!diff) {  // pattern found in substrate at position is
        Rs[0].clear() ; Rs[1].clear() ;
//        if (is==x->no_groups) err("x") ;
        if (gr_s[0]<0) { int n=-gr_s[0]-1 ; for (int i=0;i<=is;i++) Rs[n].add_right(x->groups[i]) ; }
//        cout <<"X" ;
        if (gr_s.size()>1 && gr_s[gr_s.size()-1]<0) { int n=-gr_s[gr_s.size()-1]-1 ; for (int i=x->no_groups-1;i>=is+gr_s.size()-1;i--) Rs[n].add_right(x->groups[i]) ; }
      //  if (gr_s[0]>0 && gr_s[gr_s.size()-1]>0) err("ww") ;
//        cout <<"*" ;
      
        diff=1 ;
        for (int s2=0;s2<2;s2++) {  // search for patter in the product
          int n=gr_p.size() ;
          for (int ip=0;ip<=y->no_groups-n;ip++) {
            diff=0 ;
            for (int j=0;j<n;j++) {
              if (gr_p[j]>=0 && (gr[gr_p[j]]!=y->groups[ip+j] || (grat_p[j]==1 && ip+j>0 && ip+j<y->no_groups-1)/* not at the end*/)) { diff=1 ; break ; }
            }
            if (!diff) { // pattern found in product at position ip
              // ok, basic pattern agrees, now we check if R1,R2 etc agree on both sides of equation
              ok=1 ;

 if (v)  cout <<endl<<"reaction: "<<x->string<<" --> "<<y->string<<" : " ;

//              cout <<"#"<<is<<" --> "<<ip<<"#" ;
              
              // extract the Rs
              Rp[0].clear() ; Rp[1].clear() ;
//              cout <<"$"<<gr_p[0] ;
              if (gr_p[0]<0) { int n=-gr_p[0]-1 ; for (int i=0;i<=ip;i++) Rp[n].add_right(y->groups[i]) ; }
//              cout <<"$" ;
              if (gr_p.size()>1 && gr_p[gr_p.size()-1]<0) { int n=-gr_p[gr_p.size()-1]-1 ; for (int i=y->no_groups-1;i>=ip+gr_p.size()-1;i--) Rp[n].add_right(y->groups[i]) ; }
//              cout <<"$"<<endl ;
      
  if (v) {
        cout <<"[]" ;
        for (int n=0;n<gr_s.size();n++) cout <<gr_s[n]<<" " ; cout <<endl ;
        for (int n=0;n<gr_p.size();n++) cout <<gr_p[n]<<" " ; cout <<endl ;
        for (int n=0;n<2;n++) {
      //    rs[n].to_string() ; rp[n].to_string() ; 
      //    cout <<"sub: "<<rs[n].string<<" prod: "<<rp[n].string<<endl ;
          cout <<"sub: "<<Rs[n].no_groups<<" prod: "<<Rp[n].no_groups<<endl ;
        }
      }
      
      
              Compound *Rscopy[2], *Rpcopy[2] ;
              // modify the Rs (additional/mising H's etc.)
              for (int n=0;n<2;n++) { // go through R1,R2
                Rscopy[n]=Rs[n].clone() ; Rpcopy[n]=Rp[n].clone() ; 
                for (int m=0;m<gr_s.size();m+=(gr_s.size()>1?gr_s.size()-1:1)) {
                  for (int l=0;l<gr_p.size();l+=(gr_p.size()>1?gr_p.size()-1:1)) {
//                    cout <<n<<" "<<m<<" "<<l<<endl ;
//                    cout <<gr_s[m]<<"|"<<gr_p[l]<<"|"<<string_s[m]<<"|"<<string_p[l]<<" " ;
                    if (gr_s[m]==-n-1 && gr_p[l]==-n-1 && string_s[m]=="R1" && string_p[l]=="R1H") {
                      Rscopy[n]=Rs[n].clone() ; 
if (v) {                     Rscopy[n]->to_string() ; cout <<endl<<"1!!!"<<Rscopy[n]->string<<"------>" ; }
                      Rscopy[n]->add_hydrogen_to_last() ;
if (v)  {                    Rscopy[n]->to_string() ; cout <<Rscopy[n]->string<<"!!!!!"<<endl ; }
                      Rscopy[n]->optimize(0) ;
      //          err("x") ;
                    }
                    if (gr_s[m]==-n-1 && gr_p[l]==-n-1 && string_s[m]=="R1H" && string_p[l]=="R1") {
                      Rpcopy[n]=Rp[n].clone() ; 
if (v)  {                    Rscopy[n]->to_string() ; cout <<endl<<"2!!!"<<Rscopy[n]->string<<"------>" ; }
                      Rpcopy[n]->add_hydrogen_to_last() ;
if (v)  {                    Rscopy[n]->to_string() ; cout <<Rscopy[n]->string<<"!!!!!"<<endl ;  }
                      Rpcopy[n]->optimize(0) ;
      //          err("x") ;
                    }
/*                    if (gr_s[m]==-n-1 && gr_p[l]==-n-1 && string_s[m]=="R1" && string_p[l]=="R1p") {
                      Rscopy[n]=Rs[n].clone() ; 
//                      Rscopy[n]->to_string() ; cout <<endl<<"!!!!"<<Rscopy[n]->string<<"------>" ;
                      Rscopy[n]->add_hydrogen_to_last() ;
//                      Rscopy[n]->to_string() ; cout <<Rscopy[n]->string<<"!!!!!"<<endl ;
                      Rscopy[n]->optimize(0) ;
      //          err("x") ;
                    }*/
                  }
                }
              }
              
              
              // compare the Rs 
              for (int n=0;n<2;n++) {
                if (/*rs[n].no_groups>0 && rp[n].no_groups>0 && */!equalQ(Rscopy[n],Rpcopy[n])) ok=0 ; // the Rs do not agree
              }
              
      /*        Compound *ynew=x->clone() ; //, *ycopy=y->clone() ;
              for (k=0;k<gr_s.size();k++) if (gr_s[k]>=0) ynew->groups[k+is]=gr[gr_p[k+ip]] ;
              if (gr_s[0]<0 && gr_p[0]<0) {// R1-xxx
      //          cout <<is<<endl ;
                if (string_s[0]=="R1" && string_p[0]=="R1H" && y->groups[ip]==add_hydrogen(x->groups[is],0)) ok=1 ;
              }*/
              if (ok) { // everything agrees!
                if (v) {
                  for (int n=0;n<gr_s.size();n++) cout <<gr_s[n]<<" " ; cout <<endl ;
                  for (int n=0;n<gr_p.size();n++) cout <<gr_p[n]<<" " ; cout <<endl ;
                
                  for (int n=0;n<2;n++) {
                    Rs[n].to_string() ; Rp[n].to_string() ; 
                    cout <<"sub: "<<Rs[n].string<<" prod: "<<Rp[n].string<<endl ;
                  }              
                }
                return 1 ;              
              }
                
/*              cout <<"rejected\n" ;
              for (int n=0;n<2;n++) {
              Rs[n].to_string() ; Rp[n].to_string() ; 
              cout <<"sub: "<<Rs[n].string<<" prod: "<<Rp[n].string<<endl ;
              return 0 ;
              
              }*/
        
            }
          }
          reverse(&gr_p) ; reverse(&grat_p) ; reverse(&string_p) ; // reverse order and check again
        }
      } 
    }
    reverse(&gr_s) ; reverse(&grat_s) ; reverse(&string_s) ; // reverse order and check again
  }
//  cout <<"£" ;
  return 0 ; // pattern not found in x

}

vector<Reaction_classes*> rcs ;

bool contains_xyz(string s)
{
  for (int i=0;i<s.length()-1;i++) if ((s[i]=='x' && s[i+1]=='x') || (s[i]=='y' && s[i+1]=='y') || (s[i]=='z' && s[i+1]=='z')) return true ;
  return false ; 
}

string group_match_pattern(Group *g,string p) // returns string that is a sum of all remaining groups beyond the pattern
{
//  cout <<p<<" " ;
//  cout <<g->name<<" " ;
  vector <string> ligs ;
  parse(p,&ligs,"(") ;
  if (ligs.size()==0) { // if only "C"
    if (p!="C") { cout <<p ; err("p!=Cxx && p!=Cyy") ; }
    ligs.push_back(p) ; 
  }
//  for (int i=0;i<ligs.size();i++) cout <<ligs[i]<<" " ;
//  cout <<endl ;
  if (ligs[0]!="C") err("ligs[0]!=C") ;
  if (ligs.size()-1<=g->ligands.size()) {
    int found=ligs.size()-1 ;
//    cout <<"*"<<ligs.size()<<" " ;
//    cout <<g->name<<" "<<p<<endl ; err("x") ;
    vector <bool> vis(g->ligands.size(),false) ;
    for (int i=1;i<ligs.size();i++) {
      for (int j=0;j<g->ligands.size();j++) if (vis[j]==false && ligs[i].substr(0,ligs[i].length()-1)==g->ligands[j]) { found-- ; vis[j]=true ; break ; }
    }

//    for (int i=1;i<ligs.size();i++) if (ligs[i].substr(0,ligs[i].length()-1)!=g->ligands[i-1]) ok=0 ;
    if (found==0) {
      string ns ;
      for (int j=0;j<g->ligands.size();j++) if (vis[j]==false) ns+=g->ligands[j] ;
//      cout <<g->name<<" "<<p<<endl ; 
//      cout <<ns<<endl ;
//      err("X") ;
      return ns ;
    }
  }
  return string("NULL") ; 
}

void replace_every_occurrence(string &s,string src,string dest) 
{
//  cout<<s<<": "<<src<<"----------->"<<dest<<endl ;
  int i=0;
  bool rep=false ;
  while ((i=s.find(src,i))!= string::npos) {
//    cout <<s<<endl ;
    s.replace(i,src.length(),dest) ; rep=true ;
  }
  if (!rep) err("no replacement found") ;
}

//add automatic replacement R1- and -R2 by H

void load_reaction_types() 
{
  int i ;
  string s ;
  ifstream data("reaction_types_v3_2_2.txt") ;
  while (!data.eof()) {
    getline(data,s) ; 
    if (s.length()>3) { // if more than and empty line...
      if (s[0]=='*') rcs.push_back(new Reaction_classes(s.substr(1,100))) ; // ...create new reaction class
      else if (s[0]!='#') { // or add new reaction pattern
        if (contains_xyz(s)) { // or if it contains any wildcard chars create all possible patterns
          vector <string> pr, grs, grp ; // groups to replace
          parse(s,&pr,"-=# ") ; 
          int mode=0 ; 
          for (int l=0;l<pr.size();l++) {
            if (contains_xyz(pr[l])) {
              if (mode==0) grs.push_back(pr[l].substr(0,pr[l].size()-2)) ; else grp.push_back(pr[l].substr(0,pr[l].size()-2)) ;
            }
            if (pr[l]==">") mode=1 ; // switch to products
          }
          if (grs.size()!=grp.size()) err("non-equal number of wildcard groups on s/p side") ;
          for (int l=0;l<grs.size();l++) cout <<"___"<<grs[l]<<" "<<grp[l]<<"___"<<endl ;
//          err("xxx") ;
//          groups <-- s[] identify unique groups
          
          if (grs.size()==1) {
            for (int k=0;k<_no_groups;k++) {
              string rest ;
              if (( rest=group_match_pattern(gr[k],grs[0]) )!="NULL") {
                for (int l=0;l<_no_groups;l++) {
                  if (group_match_pattern(gr[l],grp[0])==rest) {
                    string sn=s ; 
                    cout <<gr[k]->name<<"-->"<<gr[l]->name<<endl ;
                    replace_every_occurrence(sn,grs[0]+"xx",gr[k]->name) ; 
                    replace_every_occurrence(sn,grp[0]+"xx",gr[l]->name) ; 
                    cout <<endl<<sn<<endl ;
                    rcs[rcs.size()-1]->add(sn.data()) ;
                  }
                }
              }
            }
          } else if (grs.size()==2) { 
            for (int k=0;k<_no_groups;k++) {
//              cout <<grs[j].substr(0,grs[j].length()-2)<<endl ;
//              if (string(gr[k]->name).find( grs[j].substr(0,grs[j].length()-2) )!=string::npos) {
              string rest0 ;
//              cout <<"*" ;
              if (( rest0=group_match_pattern(gr[k],grs[0]) )!="NULL") {
//                if (rest=="") err("q") ;
//                cout <<gr[k]->name<<endl ;
                for (int l=0;l<_no_groups;l++) {
//                  cout <<"=" ;
//                  if (string(gr[l]->name).find( grp[j].substr(0,grp[j].length()-2) )!=string::npos) {
                  if (group_match_pattern(gr[l],grp[0])==rest0) {
//                    cout <<gr[l]->name<<endl ;
//                    if (string(gr[k]->name).compare(grs[j].length()-2,100,string(gr[l]->name),grp[j].length()-2,100)==0) { // finally, are the xx's the same in both groups?
                    string sn0=s ; 
                    cout <<sn0<<endl ;
                    cout <<gr[k]->name<<"-0->"<<gr[l]->name<<endl ;
                    replace_every_occurrence(sn0,grs[0]+"xx",gr[k]->name) ; 
                    cout <<sn0<<endl ;
                    replace_every_occurrence(sn0,grp[0]+"xx",gr[l]->name) ; 
                    cout <<sn0<<endl ;
//                    if (gr[k]->name==string("CH2") && gr[l]->name==string("CH")) err("x") ;

                    for (int m=0;m<_no_groups;m++) {
                      string rest1 ;
                      if (( rest1=group_match_pattern(gr[m],grs[1]) )!="NULL") {
                        for (int n=0;n<_no_groups;n++) {
                          if (group_match_pattern(gr[n],grp[1])==rest1) {
                            string sn1=sn0 ;
                            cout <<sn1<<endl ;
                            cout <<gr[m]->name<<"-1->"<<gr[n]->name<<endl ;
                            replace_every_occurrence(sn1,grs[1]+"yy",gr[m]->name) ; 
                            cout <<sn1<<endl ;
                            replace_every_occurrence(sn1,grp[1]+"yy",gr[n]->name) ; 
                            cout <<sn1<<endl ;
                            rcs[rcs.size()-1]->add(sn1.data()) ;
                          }
                        }
                      }
                    }
                  }          
                }
              }
            }
          } else err("too many/too few groups") ;
//          else err("no replacement rules found") ;

//          err("x") ;
        } else rcs[rcs.size()-1]->add(s.data()) ; // simply add the pattern 
      }
    }
  }
  data.close() ;    
//  err("x");
}


/*new Reaction_classes("dehydration","R1-CxH-Cx(OH)-R2  ---> R1-Cx=Cx-R2 + H2O")*/
/*Reaction_classes *rcs[3]={new Reaction_classes("tautomerism", "R1=C(OH)-R2 > R1H-CO-R2","R1=CH(OH) > R1H-CHO"),
                          new Reaction_classes("isomerisation", "R1-CO-CH2(OH) > R1-CH(OH)-CHO",  "R1-CH(OH)-CH2p > R1-CHp-CH2(OH)",
                            "R1=C(OH)-CH2p > R1=Cp-CH2(OH)", "R1-CHp-COOH > R1-CH(OH)-COp","R1=Cp-COOH > R1=C(OH)-COp"),
                  
                          new Reaction_classes("Phydrolysis", "R1-CH2p  + H2O > R1-CH2(OH) + Pi","R1-COp  + H2O > R1-COOH + Pi", "R1=CHp  + H2O > R1=CH(OH) + Pi")
                         } ;
*/


// ------------------------
// this procedure is executed from main() and has very specific arguments
void add_reaction(int type, Reaction *r, ofstream *pajek, ofstream *data, int i, int j, int &ll) 
{
  if (type) {
    //cout <<cc[i]->string<<"  --->  "<<cc[j]->string<<endl ;
    if (r==NULL) err("r==NULL") ;
// BW    cout <<r->string<<endl<<flush ;
    
    //if (j>i) 
    (*pajek)<<i+1<<" "<<j+1<<" 2 c " ; 
    if (type==1) { *pajek <<"Red"<<endl ; }
    if (type==2) { *pajek <<"Green"<<endl ; }
    if (type==3) { *pajek <<"Blue"<<endl ; }
    if (type==4) { *pajek <<"Brown"<<endl ; }
    if (type==5) { *pajek <<"Yellow"<<endl ; }        
    if (type==6) { *pajek <<"Magenta"<<endl ; }
    if (type==7) { *pajek <<"Cyan"<<endl ; }        
    if (type==8) { *pajek <<"Pink"<<endl ; }        
    if (type==9) { *pajek <<"Violet"<<endl ; }        
    if (type==10) { *pajek <<"Purple"<<endl ; }        
    if (type==11) { *pajek <<"Black"<<endl ; }        
//        data<<"\t"<<cc[j]->dG-cc[i]->dG<<endl ;
    
    int k,atp=0, nadh=0 ;
    for (k=0;k<r->no_prod;k++) if (r->prod[k]==_atp) atp++ ;
    for (k=0;k<r->no_substr;k++) if (r->substr[k]==_atp) atp-- ;
    for (k=0;k<r->no_prod;k++) if (r->prod[k]==_nadh) nadh++ ;
    for (k=0;k<r->no_substr;k++) if (r->substr[k]==_nadh) nadh-- ;

//    cout <<"$"<<r->dG<<" " ;

    *data <<r->name<<"\t"<<i<<"\t"<<j<<"\t"<<atp<<"\t"<<nadh<<"\t"<<r->dG<<"\t"<<r->string<<endl ;

    ll++ ; //cout <<"this was no."<<ll<<endl ; if (ll>=870) err("x") ; 
    delete r ;
  } else err("at this stage type should be !=0") ;
//  type=0 ;
}


//=================================================================

int main()
{
//  Compound ethanol("CH3-CH2-OH") ;
//  Compound glucose("CHO-CH(OH)-CH(OH)-CH(OH)-CH(OH)-CH2-OH") ;

//  ethanol.print() ;
//  glucose.print() ;

  int cmax=4 ;
  int p[10] ;
  int i,j,k;
  int nmax=1 ;
  for (i=0;i<cmax;i++) { p[i]=0 ; nmax*= _no_groups ; }
  
  for (i=0;i<_no_groups;i++) {
    cout <<gr[i]->name ; //<<" "<<gr[i]->bonds<<endl ; //<<gr[i]->nC<<" "<<gr[i]->nH<<" "<<gr[i]->nO<<endl ;
    if (gr[i]->bonds==1) cout <<"- " ;
    if (gr[i]->bonds==2) cout <<"= " ;
    if (gr[i]->bonds==3) cout <<"# " ;
  }
  cout <<endl<<endl ;
//  system("pause");

  save_groups() ;
  load_E_contrib() ;
  load_reaction_types() ;

  ofstream rcsfile("reaction_patterns_322.txt") ;
  for (int nr=0;nr<rcs.size();nr++) {  // reaction classes
    rcs[nr]->remove_repeated_reactions() ;
    rcs[nr]->remove_invalid_reactions() ;
    rcs[nr]->add_more_by_replacing_R_by_H() ;
    rcs[nr]->remove_repeated_reactions() ;
    for (int k=0;k<rcs[nr]->rc.size();k++) { // reaction patters within class            
      if (rcs[nr]->rc[k]->human_readable==string("")) { rcs[nr]->rc.erase(rcs[nr]->rc.begin()+k) ; k-- ; continue ; }  // BW - is this necessary????
      rcs[nr]->rc[k]->update_human_readable() ;
      rcsfile<<rcs[nr]->name<<"\t"<<rcs[nr]->rc[k]->human_readable<<endl ;
    }
  }
  rcsfile.close() ; 


  cout <<"begin:\n" ; 
  
  Compound *x=NULL ;
  Group *excl=find_group("CH(NH2)") ; // exclude this group when at the end (i.e. forming a double bond)
  int nn=0, ll=0 ;
  vector <Compound *> cc(0) ;
  for (int n=0;n<nmax;n++) { 
//    if (n%100==0) { for (j=0;j<cmax;j++) cout <<p[j]<<" " ; cout <<endl ; }
    if (1) {
    //if (gr[p[0]]->bonds==1 && gr[p[cmax-1]]->bonds==1) {
//      char xx[256],yy[256] ;
      x=new Compound() ;
      for (i=0;i<cmax;i++) x->add_right(gr[p[i]]) ;
//      x->to_string() ;
//      cout <<x->string<<"\t" ;

      if (x->dangling_bonds()==0 && x->atoms.nC>0 && (x->atoms.nO>0 || x->atoms.nP>0) && x->atoms.nH>0) {
        x->optimize() ;
        x->to_string() ;
        
        if (!find_in(x->string,"#") && x->groups[0]!=excl && x->groups[x->no_groups-1]!=excl) {// && !find_in(x->string,"=")) {   
            // reject compounds with triple bonds and "CH(NH2)=" at the end
  //        cout <<x->string<<endl ;
          if (cc.empty()) cc.push_back(x->clone()) ;
          else {
            for (i=0;i<cc.size();i++) if (equalQ(cc[i],x)) break ;
            if (i==cc.size()) cc.push_back(x->clone()) ;
          }
        }
        if (n%100==0) cout <<cc.size()<<" " ;
      } //else cout <<" rejected"<<endl ;

      //x->print() ;
    } 

    if (x!=NULL) { delete x ; x=NULL ; }

    for (j=0;j<cmax;j++) {
      p[j]++ ; if (p[j]<_no_groups) break ;
      p[j]=0 ;
    }
  }
  nn=cc.size() ;

  ofstream data, pajek("pajek_" NUM ".net") ;
  pajek<<"*Vertices "<<nn<<endl ;
    
  cout <<"\n\nlist of "<<nn<<" compounds\n" ;
  data.open("compounds_list_" NUM ".dat") ;
  for (i=0;i<nn;i++) {
    cc[i]->free_energy() ;
    cc[i]->add_name() ;
    cc[i]->if_possible_replace_dG_by_Alberty();
//    cout <<cc[i]->string<<"\t"<<cc[i]->dG<<endl ;
    data <<i<<"\t"<<cc[i]->dG<<"\t"<<cc[i]->string<<"\t" ;
    if (cc[i]->name[0]!=0) data<<cc[i]->name ; else data<<"---" ;
    data <<"\t"<<(cc[i]->charge_in_solution()) ;
    data <<endl ;
    pajek<<i+1<<" \""<<cc[i]->string<<"\\n"<<cc[i]->name<<"\""<<endl ;
  }
  data.close() ;

//  exit(0) ;
//  system("pause") ;
  
/*  cout <<"\n\nisomers\n" ;
  for (i=0;i<nn;i++) 
    for (j=i+1;j<nn;j++) 
      if (isomersQ(cc[i],cc[j])) cout <<cc[i]->string<<" and "<<cc[j]->string<<endl ;*/
      
//  dehydration(ethanol) ;


  cout <<"\n\nreactions:\n" ;
/*  for (i=0;i<nn;i++) {
    dehydration(cc[i]) ; 
  }*/
  
  pajek<<"*Edges\n" ;
  
  data.open("reactions_" NUM ".dat") ;
//  data.precision(2) ;
//  int type ;
  Reaction *r ;

  for (i=0;i<nn;i++) {
    for (j=0;j<nn;j++) {
      if (i==j) continue ;
      for (int nr=0;nr<rcs.size();nr++) {  // reaction classes
        for (int k=0;k<rcs[nr]->rc.size();k++) { // reaction patters within class
          if (rcs[nr]->rc[k]->check_if_reaction_possible(cc[i],cc[j])) {
            r=new Reaction(rcs[nr]->name.data()) ;  ;
            r->substrate(cc[i]) ; for (int l=0;l<rcs[nr]->rc[k]->ext_s.size();l++) r->substrate(rcs[nr]->rc[k]->ext_s[l]) ; 
            r->product(cc[j]) ; for (int l=0;l<rcs[nr]->rc[k]->ext_p.size();l++) r->product(rcs[nr]->rc[k]->ext_p[l]) ;
            r->update() ; 
            add_reaction(1+nr, r, &pajek, &data, i,j, ll) ;
          }
        }
      }
    }
  }
  data.close() ;
  pajek.close() ;
    
  cout <<"\n N="<<nn<<" L="<<ll<<endl ;  
  //system("pause") ;
}
