// System includes
#include <vector>
#include <iostream>

// ROOT includes
#include <TTree.h>
#include <TH1F.h>
#ifndef BRANCHWRAPPER
#define BRANCHWRAPPER

// Inherited
class BranchWrapperBase{
  public:
    BranchWrapperBase(std::string) ;
//    ~BranchWrapperBase() ;
    virtual void beginEvent() ;
    virtual void endEvent() ;
    std::string name(){ return name_; } ;
    bool  is_filled(){ return is_filled_ ; } ;
    bool is_touched(){ return is_touched_; } ;
    void touch(){ is_touched_ = true ; } ;
    void fill(){ is_filled_ = true ; touch() ; } ;
    void unfill(){ is_filled_ = false ; } ;
    virtual int config(TTree*){return -1; } ;
    virtual ~BranchWrapperBase() ;
  private:
    std::string name_ ;
    bool is_filled_ ;
    bool is_touched_ ;
};

class BranchWrapperHist  : public BranchWrapperBase{
  private:
    TH1F value_ ;
  public:
    BranchWrapperHist(std::string) ;
    ~BranchWrapperHist(){} ;
    void set(TH1F) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperB  : public BranchWrapperBase{
  private:
    bool value_ ;
  public:
    BranchWrapperB(std::string) ;
    ~BranchWrapperB(){} ;
    void set(bool) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperD  : public BranchWrapperBase{
  private:
    double value_ ;
  public:
    BranchWrapperD(std::string) ;
    ~BranchWrapperD(){} ;
    void set(double) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperF  : public BranchWrapperBase{
  private:
    float value_ ;
  public:
    BranchWrapperF(std::string) ;
    ~BranchWrapperF(){} ;
    void set(float) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperI  : public BranchWrapperBase{
  private:
    int value_ ;
  public:
    BranchWrapperI(std::string) ;
    ~BranchWrapperI(){} ;
    void set(int) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperC  : public BranchWrapperBase{
  private:
    std::string value_ ;
  public:
    BranchWrapperC(std::string) ;
    ~BranchWrapperC(){} ;
    void set(std::string) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};


class BranchWrapperU  : public BranchWrapperBase{
  private:
    unsigned int value_ ;
  public:
    BranchWrapperU(std::string) ;
    ~BranchWrapperU(){} ;
    void set(unsigned int) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperUL  : public BranchWrapperBase{
  private:
    unsigned long int value_ ;
  public:
    BranchWrapperUL(std::string) ;
    ~BranchWrapperUL(){} ;
    void set(unsigned long int) ;
    int  config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};


class BranchWrapperBV : public BranchWrapperBase{
  private:
    std::vector<bool> values_;
  public:
    BranchWrapperBV(std::string) ;
    ~BranchWrapperBV() ;
    void push(bool) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperDV : public BranchWrapperBase{
  private:
    std::vector<double> values_;
  public:
    BranchWrapperDV(std::string) ;
    ~BranchWrapperDV() ;
    void push(double) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperFV : public BranchWrapperBase{
  private:
    std::vector<float> values_;
  public:
    BranchWrapperFV(std::string) ;
    ~BranchWrapperFV() ;
    void push(float) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperIV : public BranchWrapperBase{
  private:
    std::vector<int> values_;
  public:
    BranchWrapperIV(std::string) ;
    ~BranchWrapperIV() ;
    void push(int) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperCV : public BranchWrapperBase{
  private:
    std::vector<std::string> values_;
  public:
    BranchWrapperCV(std::string) ;
    ~BranchWrapperCV() ;
    void push(std::string) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperUV : public BranchWrapperBase{
  private:
    std::vector<unsigned int> values_;
  public:
    BranchWrapperUV(std::string) ;
    ~BranchWrapperUV() ;
    void push(unsigned int) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperULV : public BranchWrapperBase{
  private:
    std::vector<unsigned long int> values_;
  public:
    BranchWrapperULV(std::string) ;
    ~BranchWrapperULV() ;
    void push(unsigned long int) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};


class BranchWrapperBVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<bool> > values_;
  public:
    BranchWrapperBVV(std::string) ;
    ~BranchWrapperBVV() ;
    void push(std::vector<bool>) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperDVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<double> > values_;
  public:
    BranchWrapperDVV(std::string) ;
    ~BranchWrapperDVV() ;
    void push(std::vector<double>) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperFVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<float> > values_;
  public:
    BranchWrapperFVV(std::string) ;
    ~BranchWrapperFVV() ;
    void push(std::vector<float>) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperIVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<int> > values_;
  public:
    BranchWrapperIVV(std::string) ;
    ~BranchWrapperIVV() ;
    void push(std::vector<int>) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};

class BranchWrapperUVV: public BranchWrapperBase{
  private:
    std::vector<std::vector<unsigned int> > values_;
  public:
    BranchWrapperUVV(std::string) ;
    ~BranchWrapperUVV() ;
    void push(std::vector<unsigned int>) ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};


// Templated (not used, yet)
template <class T>
class branchWrapper_simple{
    std::string name_ ;
    T value_ ;
  public:
    std::string name(){ return name_ ; } ;
    branchWrapper_simple(std::string) ;
    ~branchWrapper_simple(){} ;
    void set(T) ;
    T get() ;
    int config(TTree*) ;
    void beginEvent() ;
    void endEvent() ;
};
template <class T>
class branchWrapper_vector{
    std::string name_ ;
    std::vector<T> values_ ;
  public:
    std::string name(){ return name_ ; } ;
    branchWrapper_vector(std::string) ;
    ~branchWrapper_vector(){} ;
    bool push(T) ;
    T get(int) ;
    int nEntries() ;
    int config(TTree*) ;
    bool beginEvent() ;
    bool endEvent() ;
};


#endif
