using Proteomics.PSM;
using Readers;
using Readers.QuantificationResults;

namespace Plots;

public class Plot
{

}

public class Results
{
    public string ResultsPath { get; set; }
    public List<QuantifiedPeak> AllQuantifiedPeaks { get; set; }
    public List<PsmFromTsv> AllPsms { get; set; }
    public List<PsmFromTsv> AllPeptides { get; set; }
    public Dictionary<string, (PsmFromTsvFile Peptides, PsmFromTsvFile PSMs, PsmFromTsvFile QuantifiedPeaks)> IndividualFileResults { get; set; }

    public Results(string searchTaskPath)
    {
        ResultsPath = searchTaskPath;

        var allQuantifiedPeaksFile = new QuantifiedPeakFile(Path.Join(searchTaskPath, "AllQuantifiesPeaks.tsv"));
        allQuantifiedPeaksFile.LoadResults();
        AllQuantifiedPeaks = allQuantifiedPeaksFile.Results;
        
        var allPsmsFile = new PsmFromTsvFile(Path.Join(searchTaskPath, "AllPSMs.psmtsv"));
        allPsmsFile.LoadResults();
        AllPsms = allPsmsFile.Results;

        var allPeptidesFile = new PsmFromTsvFile(Path.Join(searchTaskPath, "AllPeptides.psmtsv"));
        allPeptidesFile.LoadResults();
        AllPeptides = allPeptidesFile.Results;

        var filesUsed = AllPsms.Select(x => x.FileNameWithoutExtension)
            .Distinct()
            .ToList();

        var filesInIndividualResults = Directory.GetFiles(Path.Join(searchTaskPath, "Individual File Results"));

        foreach (var file in filesUsed)
        {

        }
    }
}