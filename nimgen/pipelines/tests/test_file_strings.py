"""Provide tests for all constant strings to generate files."""


def test_file_strings():
    """Test the string constants needed to generate files for pipeline."""
    from nimgen.pipelines import (
        _htcondor_python_strings,
        _htcondor_submit_strings,
    )

    assert isinstance(_htcondor_submit_strings.TEMPLATE_DAG, str)
    assert isinstance(
        _htcondor_submit_strings.TEMPLATE_JOB.format(
            "/path/to/initial_dir", 10, 20, 30
        ),
        str,
    )
    assert isinstance(
        _htcondor_submit_strings.TEMPLATE_QUEUED_JOB.format(
            "./run_in_venv.sh step3 my_args", logs="/path/to/logs", job_id="5"
        ),
        str,
    )
    assert isinstance(
        _htcondor_python_strings.RUN_IN_VENV.format(venv="venv"), str
    )
