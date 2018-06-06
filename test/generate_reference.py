import os

from test import clean_output, run_test

if __name__ == "__main__":

    output_dir = "output"
    clean_output(output_dir)

    tests = [f for f in os.listdir("../demo/") if f.endswith(".py")]
    print("Found %d tests" % len(tests))

    # Step into output directory
    os.chdir(output_dir)

    failures = []

    for test in tests:
        run_test(test, generate_reference=1)
