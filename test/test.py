__author__ = "Marie E. Rognes (meg@simula.no)"
__copyright__ = "Copyright (C) 2010 -- Marie E. Rognes"
__license__  = "GNU GPL version 3 or any later version"

# Note: Largely copied from ffc tests

import os, re, sys, shutil, subprocess, difflib

# Global log file
logfile = None

def clean_output(dir):
    "Clean out old output directory"
    if os.path.isdir(dir):
        shutil.rmtree(dir)
    os.mkdir(dir)

def run_command(command):
    "Run command and collect errors in log file."
    (status, output) = subprocess.getstatusoutput(command)
    if status == 0:
        return (True, output)
    global logfile
    if logfile is None:
        logfile = open("../error.log", "w")
    logfile.write(output + "\n")
    return (False, output)

def run_test(test, generate_reference=False):

    test_base = test.split(".")[0]

    print("-"*80)
    print("Checking ", test_base)
    print("-"*80)

    # Copy test into this directory
    shutil.copy("../../demo/%s" % test, ".")

    # Run test
    print("Running %s..." % test, end=' ')
    (ok, output) = run_command("python %s" % test)
    if not ok:
        print("failed.")
        print(output)
        return test
    else:
        print("ok!")

    # Check that test matches old results
    reference_filename = "../references/%s.log" % test_base

    if generate_reference:
        print("Generating reference for %s" % test)
        reference_file = open(reference_filename, 'w')
        reference_file.write(output)
        return test

    if os.path.isfile(reference_filename):
        reference = open(reference_filename).read()
    else:
        print("Missing reference for %s" % test)
        return

    # Match output and reference
    if output == reference:
        print("%s-output matches reference" % test_base)
    else:
        print("%s-output differs from reference" % test_base)
        diff = "\n".join([line for line in difflib.unified_diff(reference.split("\n"), output.split("\n"))])
        print("diff = ", diff)


if __name__ == "__main__":

    output_dir = "output"
    clean_output(output_dir)

    tests = [f for f in os.listdir("../demo/") if f.endswith(".py")]
    print("Found %d tests" % len(tests))

    # Step into output directory
    os.chdir(output_dir)

    failures = []

    for test in tests:
        run_test(test)
