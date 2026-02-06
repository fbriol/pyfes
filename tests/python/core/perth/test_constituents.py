# Copyright (c) 2025 CNES
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""Comprehensive tests for the pyfes.core constituent functions (Perth)."""

from pyfes import core
import pytest


class TestConstituents:
    """Test suite for constituent parsing and querying functions (Perth)."""

    def test_parse_constituent_returns_int(self) -> None:
        """Test that parse_constituent() returns an integer ID."""
        result = core.parse_constituent('M2')
        assert isinstance(result, int)
        assert result >= 0

    def test_parse_constituent_major_waves(self) -> None:
        """Test parsing major tidal wave names."""
        major_constituents = ['M2', 'S2', 'K1', 'O1', 'N2', 'K2']

        for name in major_constituents:
            result = core.parse_constituent(name)
            assert isinstance(result, int)
            assert result >= 0

    def test_parse_constituent_invalid(self) -> None:
        """Test that parsing an invalid constituent raises an error."""
        with pytest.raises(ValueError):
            core.parse_constituent('__INVALID__')

    def test_known_constituents_returns_list(self) -> None:
        """Test that known_constituents() returns a list of strings."""
        known = core.known_constituents()

        assert isinstance(known, list)
        assert len(known) > 0
        for name in known:
            assert isinstance(name, str)
            assert len(name) > 0

    def test_known_constituents_contains_major_waves(self) -> None:
        """Test that known constituents include major tidal waves."""
        known = core.known_constituents()

        major = ['M2', 'S2', 'K1', 'O1', 'N2', 'K2']
        for name in major:
            assert name in known, f'{name} not found in known constituents'

    def test_known_constituents_size(self) -> None:
        """Test that known constituents has a reasonable size."""
        known = core.known_constituents()
        # Should have at least the major constituents
        assert len(known) >= 6

    def test_known_constituents_all_parseable(self) -> None:
        """Test that all known constituent names can be parsed."""
        known = core.known_constituents()

        for name in known:
            result = core.parse_constituent(name)
            assert isinstance(result, int)
            assert result >= 0

    def test_parse_constituent_unique_ids(self) -> None:
        """Test that different constituents have different IDs."""
        m2_id = core.parse_constituent('M2')
        s2_id = core.parse_constituent('S2')
        k1_id = core.parse_constituent('K1')

        assert m2_id != s2_id
        assert m2_id != k1_id
        assert s2_id != k1_id

    def test_parse_constituent_consistency(self) -> None:
        """Test that repeated calls return consistent results."""
        id1 = core.parse_constituent('M2')
        id2 = core.parse_constituent('M2')

        assert id1 == id2

    def test_known_constituents_consistency(self) -> None:
        """Test that multiple calls to known_constituents return same data."""
        known1 = core.known_constituents()
        known2 = core.known_constituents()

        assert len(known1) == len(known2)
        assert known1 == known2

    def test_constituents_wave_properties_via_table(self) -> None:
        """Test accessing wave properties via Perth WaveTable."""
        wt = core.perth.WaveTable(['M2', 'K1', 'S2'])

        m2 = wt['M2']

        # Should have frequency
        assert m2.frequency() > 0

        # Should have type
        assert m2.type is not None

        # Should have name
        assert m2.name == 'M2'

    def test_constituents_major_constituents_present(self) -> None:
        """Test that all major tidal constituents are present in WaveTable."""
        wt = core.perth.WaveTable()

        major = ['M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1', 'Q1']

        for name in major:
            assert name in wt, f'{name} not found in WaveTable'
