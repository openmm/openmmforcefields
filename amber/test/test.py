import os
import unittest

os.chdir("..")


class Testamber2ommmScript(unittest.TestCase):
    def test_leaprc(self):
        """Test conversion of a leaprc"""
        cmd = (
            "python amber2openmm.py -i test/leaprc.ff14SB -if leaprc -od test/ffxml "
            "--no-log --protein-test --nucleic-test"
        )
        os.system(cmd)

    def test_yaml(self):
        """Test conversion of a short yaml"""
        cmd = "python amber2openmm.py -i test/test.yaml -od test/ffxml --no-log"
        os.system(cmd)


if __name__ == "__main__":
    unittest.main()
