# -*- coding: utf-8 -*-

"""
This module provides the XMLTestRunner class, which is heavily based on the
default TextTestRunner.
"""

import os
import sys
import time
try:
    from unittest2.runner import TextTestRunner
    from unittest2.runner import TextTestResult as _TextTestResult
    from unittest2.result import TestResult
except ImportError:
    from unittest import TestResult, _TextTestResult, TextTestRunner

try:
    # Removed in Python 3
    from cStringIO import StringIO
except ImportError:
    from io import StringIO


# Allow version to be detected at runtime.
from .version import __version__, __version_info__


class _DelegateIO(object):
    """
    This class defines an object that captures whatever is written to
    a stream or file.
    """

    def __init__(self, delegate):
        self._captured = StringIO()
        self.delegate = delegate

    def write(self, text):
        self._captured.write(text)
        self.delegate.write(text)

    def __getattr__(self, attr):
        return getattr(self._captured, attr)


def testcase_name(test_method):
    testcase = type(test_method)

    # Ignore module name if it is '__main__'
    module = testcase.__module__ + '.'
    if module == '__main__.':
        module = ''
    result = module + testcase.__name__
    return result


class _TestInfo(object):
    """
    This class keeps useful information about the execution of a
    test method.
    """

    # Possible test outcomes
    (SUCCESS, FAILURE, ERROR, SKIP) = range(4)

    def __init__(self, test_result, test_method, outcome=SUCCESS, err=None):
        self.test_result = test_result
        self.test_method = test_method
        self.outcome = outcome
        self.elapsed_time = 0
        self.err = err

        self.test_description = self.test_result.getDescription(test_method)
        self.test_exception_info = (
            '' if outcome in (self.SUCCESS, self.SKIP)
            else self.test_result._exc_info_to_string(
                    self.err, test_method)
        )

        self.test_name = testcase_name(test_method)
        self.test_id = test_method.id()

    def id(self):
        return self.test_method.id()

    def test_finished(self):
        """Save info that can only be calculated once a test has run.
        """
        self.elapsed_time = \
            self.test_result.stop_time - self.test_result.start_time

    def get_description(self):
        """
        Return a text representation of the test method.
        """
        return self.test_description

    def get_error_info(self):
        """
        Return a text representation of an exception thrown by a test
        method.
        """
        return self.test_exception_info


class _XMLTestResult(_TextTestResult):
    """
    A test result class that can express test results in a XML report.

    Used by XMLTestRunner.
    """
    def __init__(self, stream=sys.stderr, descriptions=1, verbosity=1,
                 elapsed_times=True):
        _TextTestResult.__init__(self, stream, descriptions, verbosity)
        self.successes = []
        self.callback = None
        self.elapsed_times = elapsed_times

    def _prepare_callback(self, test_info, target_list, verbose_str,
                          short_str):
        """
        Appends a _TestInfo to the given target list and sets a callback
        method to be called by stopTest method.
        """
        target_list.append(test_info)

        def callback():
            """Prints the test method outcome to the stream, as well as
            the elapsed time.
            """

            test_info.test_finished()

            # Ignore the elapsed times for a more reliable unit testing
            if not self.elapsed_times:
                self.start_time = self.stop_time = 0

            if self.showAll:
                self.stream.writeln(
                    '%s (%.3fs)' % (verbose_str, test_info.elapsed_time)
                )
            elif self.dots:
                self.stream.write(short_str)
        self.callback = callback

    def startTest(self, test):
        """
        Called before execute each test method.
        """
        self.start_time = time.time()
        TestResult.startTest(self, test)

        if self.showAll:
            self.stream.write('  ' + self.getDescription(test))
            self.stream.write(" ... ")

    def stopTest(self, test):
        """
        Called after execute each test method.
        """
        _TextTestResult.stopTest(self, test)
        self.stop_time = time.time()

        if self.callback and callable(self.callback):
            self.callback()
            self.callback = None

    def addSuccess(self, test):
        """
        Called when a test executes successfully.
        """
        self._prepare_callback(
            _TestInfo(self, test), self.successes, 'OK', '.'
        )

    def addFailure(self, test, err):
        """
        Called when a test method fails.
        """
        testinfo = _TestInfo(self, test, _TestInfo.ERROR, err)
        self.errors.append((
            testinfo,
            self._exc_info_to_string(err, test)
        ))
        self._prepare_callback(testinfo, [], 'FAIL', 'F')

    def addError(self, test, err):
        """
        Called when a test method raises an error.
        """
        testinfo = _TestInfo(self, test, _TestInfo.ERROR, err)
        self.errors.append((
            testinfo,
            self._exc_info_to_string(err, test)
        ))
        self._prepare_callback(testinfo, [], 'ERROR', 'E')

    def addSkip(self, test, reason):
        """
        Called when a test method was skipped.
        """
        testinfo = _TestInfo(self, test, _TestInfo.SKIP, reason)
        self.skipped.append((testinfo, reason))
        self._prepare_callback(testinfo, [], 'SKIP', 'S')

    def printErrorList(self, flavour, errors):
        """
        Writes information about the FAIL or ERROR to the stream.
        """
        for test_info, error in errors:
            self.stream.writeln(self.separator1)
            self.stream.writeln(
                '%s [%.3fs]: %s' % (flavour, test_info.elapsed_time,
                                    test_info.get_description())
            )
            self.stream.writeln(self.separator2)
            self.stream.writeln('%s' % test_info.get_error_info())

    def _get_info_by_testcase(self, outsuffix):
        """
        Organizes test results by TestCase module. This information is
        used during the report generation, where a XML report will be created
        for each TestCase.
        """
        tests_by_testcase = {}

        for tests in (self.successes, self.failures, self.errors, self.skipped):
            for test_info in tests:
                if isinstance(test_info, tuple):
                    # This is a skipped, error or a failure test case
                    test_info = test_info[0]
                testcase_name = test_info.test_name
                if not testcase_name in tests_by_testcase:
                    tests_by_testcase[testcase_name] = []
                tests_by_testcase[testcase_name].append(test_info)

        return tests_by_testcase

    def _report_testsuite(suite_name, outsuffix, tests, xml_document):
        """
        Appends the testsuite section to the XML document.
        """
        testsuite = xml_document.createElement('testsuite')
        xml_document.appendChild(testsuite)

        testsuite.setAttribute('name', "%s-%s" % (suite_name, outsuffix))
        testsuite.setAttribute('tests', str(len(tests)))

        testsuite.setAttribute(
            'time', '%.3f' % sum(map(lambda e: e.elapsed_time, tests))
        )
        failures = filter(lambda e: e.outcome == _TestInfo.FAILURE, tests)
        testsuite.setAttribute('failures', str(len(list(failures))))

        errors = filter(lambda e: e.outcome == _TestInfo.ERROR, tests)
        testsuite.setAttribute('errors', str(len(list(errors))))

        return testsuite

    _report_testsuite = staticmethod(_report_testsuite)

    def _test_method_name(test_id):
        """
        Returns the test method name.
        """
        return test_id.split('.')[-1]

    _test_method_name = staticmethod(_test_method_name)

    def _report_testcase(suite_name, test_result, xml_testsuite, xml_document):
        """
        Appends a testcase section to the XML document.
        """
        testcase = xml_document.createElement('testcase')
        xml_testsuite.appendChild(testcase)

        testcase.setAttribute('classname', suite_name)
        testcase.setAttribute(
            'name', _XMLTestResult._test_method_name(test_result.test_id)
        )
        testcase.setAttribute('time', '%.3f' % test_result.elapsed_time)

        if (test_result.outcome != _TestInfo.SUCCESS):
            elem_name = ('failure', 'error', 'skipped')[test_result.outcome - 1]
            failure = xml_document.createElement(elem_name)
            testcase.appendChild(failure)
            if test_result.outcome != _TestInfo.SKIP:
                failure.setAttribute('type', test_result.err[0].__name__)
                failure.setAttribute('message', str(test_result.err[1]))
                error_info = str(test_result.get_error_info())
                failureText = xml_document.createCDATASection(error_info)
                failure.appendChild(failureText)
            else:
                failure.setAttribute('type', 'skip')
                failure.setAttribute('message', test_result.err)


    _report_testcase = staticmethod(_report_testcase)

    def _report_output(test_runner, xml_testsuite, xml_document):
        """
        Appends the system-out and system-err sections to the XML document.
        """
        systemout = xml_document.createElement('system-out')
        xml_testsuite.appendChild(systemout)

        systemout_text = xml_document.createCDATASection(sys.stdout.getvalue())
        systemout.appendChild(systemout_text)

        systemerr = xml_document.createElement('system-err')
        xml_testsuite.appendChild(systemerr)

        systemerr_text = xml_document.createCDATASection(sys.stderr.getvalue())
        systemerr.appendChild(systemerr_text)

    _report_output = staticmethod(_report_output)

    def generate_reports(self, test_runner):
        """
        Generates the XML reports to a given XMLTestRunner object.
        """
        from xml.dom.minidom import Document
        all_results = self._get_info_by_testcase(test_runner.outsuffix)

        if (isinstance(test_runner.output, str) and not
                os.path.exists(test_runner.output)):
            os.makedirs(test_runner.output)

        for suite, tests in all_results.items():
            doc = Document()

            # Build the XML file
            testsuite = _XMLTestResult._report_testsuite(
                suite, test_runner.outsuffix, tests, doc
            )
            for test in tests:
                _XMLTestResult._report_testcase(suite, test, testsuite, doc)
            _XMLTestResult._report_output(test_runner, testsuite, doc)
            xml_content = doc.toprettyxml(indent='\t')

            if type(test_runner.output) is str:
                report_file = open(
                    '%s%sTEST-%s-%s.xml' % (
                        test_runner.output, os.sep, suite,
                        test_runner.outsuffix
                    ), 'w'
                )
                try:
                    report_file.write(xml_content)
                finally:
                    report_file.close()
            else:
                # Assume that test_runner.output is a stream
                test_runner.output.write(xml_content)


class XMLTestRunner(TextTestRunner):
    """
    A test runner class that outputs the results in JUnit like XML files.
    """
    def __init__(self, output='.', outsuffix=None, stream=sys.stderr,
                 descriptions=True, verbosity=1, elapsed_times=True):
        TextTestRunner.__init__(self, stream, descriptions, verbosity)
        self.verbosity = verbosity
        self.output = output
        if outsuffix:
            self.outsuffix = outsuffix
        else:
            self.outsuffix = time.strftime("%Y%m%d%H%M%S")
        self.elapsed_times = elapsed_times

    def _make_result(self):
        """
        Creates a TestResult object which will be used to store
        information about the executed tests.
        """
        return _XMLTestResult(
            self.stream, self.descriptions, self.verbosity, self.elapsed_times
        )

    def _patch_standard_output(self):
        """
        Replaces stdout and stderr streams with string-based streams
        in order to capture the tests' output.
        """
        sys.stdout = _DelegateIO(sys.stdout)
        sys.stderr = _DelegateIO(sys.stderr)

    def _restore_standard_output(self):
        """
        Restores stdout and stderr streams.
        """
        sys.stdout = sys.stdout.delegate
        sys.stderr = sys.stderr.delegate

    def run(self, test):
        """
        Runs the given test case or test suite.
        """
        try:
            # Prepare the test execution
            self._patch_standard_output()
            result = self._make_result()

            # Print a nice header
            self.stream.writeln()
            self.stream.writeln('Running tests...')
            self.stream.writeln(result.separator2)

            # Execute tests
            start_time = time.time()
            test(result)
            stop_time = time.time()
            time_taken = stop_time - start_time

            # Print results
            result.printErrors()
            self.stream.writeln(result.separator2)
            run = result.testsRun
            self.stream.writeln("Ran %d test%s in %.3fs" % (
                run, run != 1 and "s" or "", time_taken)
            )
            self.stream.writeln()

            expectedFails = unexpectedSuccesses = skipped = 0
            try:
                results = map(len, (result.expectedFailures,
                                    result.unexpectedSuccesses,
                                    result.skipped))
            except AttributeError:
                pass
            else:
                expectedFails, unexpectedSuccesses, skipped = results

            # Error traces
            infos = []
            if not result.wasSuccessful():
                self.stream.write("FAILED")
                failed, errored = map(len, (result.failures, result.errors))
                if failed:
                    infos.append("failures={0}".format(failed))
                if errored:
                    infos.append("errors={0}".format(errored))
            else:
                self.stream.write("OK")

            if skipped:
                infos.append("skipped={0}".format(skipped))
            if expectedFails:
                infos.append("expected failures={0}".format(expectedFails))
            if unexpectedSuccesses:
                infos.append("unexpected successes={0}".fornat(unexpectedSuccesses))

            if infos:
                self.stream.writeln(" ({0})".format(", ".join(infos)))
            else:
                self.stream.write("\n")

            # Generate reports
            self.stream.writeln()
            self.stream.writeln('Generating XML reports...')
            result.generate_reports(self)
        finally:
            self._restore_standard_output()

        return result
