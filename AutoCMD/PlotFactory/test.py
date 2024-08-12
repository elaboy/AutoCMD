import utils 
import unittest

class TestingUtils(unittest.TestCase):
    def test_get_dataframes(self):
        paths_to_psms = [r"E:\Analyzed\A549_Full_Chronologer_2\Task4SearchTask\AllPSMs.psmtsv",
                         r"E:\Analyzed\A549_Full_Chronologer_2\Task4SearchTask\AllPSMs.psmtsv",
                         r"E:\Analyzed\A549_Full_Chronologer_2\Task4SearchTask\AllPSMs.psmtsv",
                         r"E:\Analyzed\A549_Full_Chronologer_2\Task4SearchTask\AllPSMs.psmtsv"]
        
        self.assertEqual(len(utils.get_dataframes(paths_to_psms)), 4)

if __name__ == "__main__":
    unittest.main()