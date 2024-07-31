using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Security.Cryptography;
using System.Text;
using System.Text.RegularExpressions;

namespace AutoCMD;

public class CMD
{
    public List<string> Spectra { get; set; }
    public List<string> Tasks { get; set; }
    public List<string> Databases { get; set; }
    public string? MMCmdPath { get; set; }
    public string OutputPath { get; set; }
    public string Runner { get; set; }

    public static void Main(string[] args)
    {
        var cmd = new CMD(args[0], args[1],
            int.Parse(args[2]), args[3]);

        StringBuilder taskStringBuilder = new StringBuilder();
        taskStringBuilder.Append("-t ");
        foreach (var task in cmd.Tasks)
        {
            taskStringBuilder.Append(task + " ");
        }
        string tasks = taskStringBuilder.ToString();

        StringBuilder spectraStringBuilder = new StringBuilder();
        spectraStringBuilder.Append("-s ");
        foreach (var spectra in cmd.Spectra)
        {
            spectraStringBuilder.Append(spectra + " ");
        }
        string spectras = spectraStringBuilder.ToString();


        StringBuilder databaseStringBuilder = new StringBuilder();
        databaseStringBuilder.Append("-d ");
        foreach (var database in cmd.Databases)
        {
            databaseStringBuilder.Append(database + " ");
        }

        string databases = databaseStringBuilder.ToString();

        string arguments = $"{spectras} {tasks} {databases} -o {cmd.OutputPath}";
        try
        {
            Process runMM = new Process();
            runMM.StartInfo.UseShellExecute = false;
            runMM.StartInfo.FileName = cmd.MMCmdPath;
            runMM.StartInfo.Arguments = arguments;
            runMM.StartInfo.RedirectStandardOutput = true;
            runMM.StartInfo.RedirectStandardError = true;
            runMM.StartInfo.CreateNoWindow = false;
            runMM.OutputDataReceived += (sender, e) => Console.WriteLine(e.Data);
            runMM.ErrorDataReceived += (sender, e) => Console.WriteLine("Error: " + e.Data);

            runMM.Start();
            runMM.BeginOutputReadLine();
            runMM.BeginErrorReadLine();
            runMM.WaitForExit();

            Console.WriteLine("Process completed with exit code: " + runMM.ExitCode);
        }
        catch(Exception e)
        {
            Console.WriteLine(e.Message);
        }
    }
    public CMD(string runner, string setName, int maxThreads, string outputPath)
    {
        MMCmdPath = GetMMPath(runner);
        GetFiles(setName, maxThreads);
        Databases = new List<string>()
        {
            @"E:\Proteomes\uniprotkb_proteome_UP000005640_AND_revi_2024_07_08.xml" 
        };
        OutputPath = outputPath;
    }

    static string GetMMPath(string runnerType)
    {
        switch (runnerType)
        {
            case "Experimental":
                return @"E:\MetaMorpheusVersions\MetaMorpheus\CMD\bin\Release\net6.0\CMD.exe";
            //for the moment
            default:
                return @"E:\MetaMorpheusVersions\MetaMorpheus\CMD\bin\Release\net6.0\CMD.exe"; 
        }
    }

    void GetFiles(string setName, int threadsToUse)
    {
        var spectraBaseDirectory = @"E:\Datasets\Mann_11cell_lines";

        List<string> filePaths = new();

        var path = Path.Join(spectraBaseDirectory, setName+@"\");

        foreach (var directory in Directory.GetDirectories(path))
        {
            var rawFiles = Directory.GetFiles(directory, "*.raw");
            filePaths.AddRange(rawFiles);
        }

        Spectra = filePaths;

        Tasks = new();
        // get task tomls 
        foreach (var task in Directory.GetFiles(@"E:\TaskSettings"))
        {
            // enter the file to add the desired threads to use
            StreamReader reader = new StreamReader(task);
            string input = reader.ReadToEnd();
            reader.Close();

            string output = ReplaceMaxThreads(input, threadsToUse);

            File.WriteAllText(task, output);

            Tasks.Add(task);
        }
    }

    static string ReplaceMaxThreads(string input, int threadsToUse)
    {
        string pattern = @"MaxThreadsToUsePerFile\s*=\s*\d+";
        string replacement = $"MaxThreadsToUsePerFile = {threadsToUse}";
        return Regex.Replace(input, pattern, replacement);
    }
}