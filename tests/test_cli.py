from click.testing import CliRunner
from tfbs.cli import cli
import pytest

@pytest.mark.skip(reason="no way of currently testing this")
def test_tfbs_cli():
    runner = CliRunner()
    result = runner.invoke(cli, ["--fasta", "test.fa", "--pwms", "test_pwm.json"])
    assert result.exit_code == 0