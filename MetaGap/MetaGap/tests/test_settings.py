import os
from unittest.mock import patch

from django.test import SimpleTestCase

from MetaGap.settings import _split_env_list


class SplitEnvListTests(SimpleTestCase):
    def test_split_env_list_strips_and_filters_values(self) -> None:
        with patch.dict(os.environ, {"TEST_LIST": " a , ,b ,c "}, clear=True):
            self.assertEqual(_split_env_list("TEST_LIST", ""), ["a", "b", "c"])

    def test_split_env_list_falls_back_to_default(self) -> None:
        with patch.dict(os.environ, {}, clear=True):
            self.assertEqual(_split_env_list("TEST_LIST", "foo,bar"), ["foo", "bar"])
