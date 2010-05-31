#ifndef TagProbeFitter_h
#define TagProbeFitter_h

#include "TFile.h"
#include "TChain.h"
#include "TGraphAsymmErrors.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooDataSet.h"

class TagProbeFitter: public TGraphAsymmErrors{
  public:
  ///construct the fitter with the inputFileName, inputDirectoryName, inputTreeName, outputFileName and specify wether to save the workspace with data for each bin 
  TagProbeFitter(std::vector<std::string> inputFileNames, std::string inputDirectoryName, std::string inputTreeName, std::string outputFileName, int numCPU = 1, bool saveWorkspace_ = false, bool floatShapeParameters = true, std::vector<std::string> fixVars_ = std::vector<std::string>() );

  ///destructor closes the files
  ~TagProbeFitter();

  ///adds a new real variable to the set of variables describing the data in the tree
  bool addVariable(std::string variableName, std::string title, double low, double hi, std::string units);

  ///adds a new category variable to the set of variables describing the data in the tree; "expression" is parsed by factory()
  bool addCategory(std::string categoryName, std::string title, std::string expression);

  ///add a new PDF to the list of available PDFs; "pdfCommands" are parsed by factory().
  /// the user needs to define efficiency[0.9,0,1] for the initial value, "signal" PDF, "backgroundPass" PDF and "backgroundFail" PDF
  void addPdf(std::string pdfName, std::vector<std::string>& pdfCommands);

  ///set a list of variables to fix during first fit iteration. If the list is empty, do one iteration.
  void addFixedVariavles(std::vector<string>);

  ///calculate the efficiency for a particular binning of the data; it saves everything in the directory "dirName", uses the previously defined PDF with name "pdfName"
  std::string calculateEfficiency(std::string dirName, std::string efficiencyCategory, std::string efficiencyState, std::vector<std::string>& unbinnedVariables, std::map<std::string, std::vector<double> >& binnedReals, std::map<std::string, std::vector<std::string> >& binnedCategories, std::vector<std::string>& binToPDFmap, bool saveWork);

  /// set number of bins to use when making the plots; 0 = automatic
  void setBinsForMassPlots(int bins) ;

  protected:
  ///pointer to the input TTree Chain of data
  TChain* inputTree;

  ///pointer to the output file
  TFile* outputFile;

  ///pointer to the TDirectory in the output file that is the root directory for this fitter
  TDirectory* outputDirectory;

  ///number of CPUs to use for the fit
  int numCPU;

  ///the default option wether to save the workspace for each bin
  bool saveWorkspace;

  ///number of bins to use in mass shape plots; 0 = automatic
  int massBins;

  ///the map of pdf names to the vector of commands to build the pdf
  std::map<std::string, std::vector<std::string> > pdfs;

  ///the set of variables describing the data in the input TTree
  RooArgSet variables;

  ///list of variables fo fix (see below)
  std::vector<std::string> fixVars;
  std::vector<double> fixVarValues;

  ///release some variables before the fit in each bin
  ///if set to "false" will fit all dataset to get values of specified variables and then fit all bins having them fixed
  ///if set to "true" (default) will not fit all dataset, just each bin with fixed and then released variables
  bool floatShapeParameters;

  ///a RooWorkspace object to parse input parameters with ".factory()"
  RooWorkspace parameterParser;

  ///fix or release variables selected by user
  void varFixer(RooWorkspace* w, bool fix);
  ///store values in the vector
  void varSaver(RooWorkspace* w);
  ///restore variables's values for fit starting point
  void varRestorer(RooWorkspace* w);

  ///calculate the efficiecny with a simulataneous maximum likelihood fit in the dataset found in the workspace with PDF pdfName
  void doFitEfficiency(RooWorkspace* w, std::string pdfName, RooRealVar& efficiency);

  ///calculate the efficiecny with side band substraction in the dataset found in the workspace
  void doSBSEfficiency(RooWorkspace* w, RooRealVar& efficiency);

  ///calculate the efficiecny by counting in the dataset found in the workspace
  void doCntEfficiency(RooWorkspace* w, RooRealVar& efficiency);

  ///creates the simultaneous PDF in the workspace according to the "pdfCommands"
  void createPdf(RooWorkspace* w, std::vector<std::string>& pdfCommands);

  ///sets initial values of the PDF parameters based on the data available in the workspace
  void setInitialValues(RooWorkspace* w);

  ///saves the fit canvas
  void saveFitPlot(RooWorkspace* w);

  ///saves the distributions canvas
  void saveDistributionsPlot(RooWorkspace* w);

  ///saves the efficiency plots
  void saveEfficiencyPlots(RooDataSet& eff, TString effName, RooArgSet& binnedVariables, RooArgSet& mappedCategories);
  
  ///makes the 1D plot
  void makeEfficiencyPlot1D(RooDataSet& eff, RooRealVar& v, TString plotName, TString plotTitle, TString effName);
  
  ///makes the 2D plot
  void makeEfficiencyPlot2D(RooDataSet& eff, RooRealVar& v1, RooRealVar& v2, TString plotName, TString plotTitle, TString effName);
  
};

#endif //TagProbeFitter_h
