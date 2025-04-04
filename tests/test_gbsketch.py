"""
Tests for gbsketch
"""
import os
import csv
import pytest

import sourmash
from sourmash import sourmash_args
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed



def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch')

    assert 'usage:  gbsketch' in runtmp.last_result.err


def test_gbsketch_simple(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    print(f"looking for path: {output}")

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()
    assert os.path.exists(failed)
    with open(failed, 'r') as failF:
        fail_lines = failF.readlines()
        print(fail_lines)
        assert len(fail_lines) == 2
        assert fail_lines[0] == "accession,name,moltype,md5sum,download_filename,url,range\n"
        acc, name, moltype, md5sum, download_filename, url, range = fail_lines[1].strip().split(',')
        assert acc == "GCA_000175535.1"
        assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
        assert moltype == "protein"
        assert download_filename == "GCA_000175535.1_protein.faa.gz"
        assert url == ""
        assert range == ""


def test_gbsketch_simple_default_failed(runtmp, capfd):
    # test the default value for --failed
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('acc.csv.fail.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '-r', '3',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200",
                    in_dir=runtmp.output(''))

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    print(f"looking for path: {output}")

    assert os.path.exists(failed)
    with open(failed, 'r') as failF:
        fail_lines = failF.readlines()
        print(fail_lines)
        assert len(fail_lines) == 2
        assert fail_lines[0] == "accession,name,moltype,md5sum,download_filename,url,range\n"
        acc, name, moltype, md5sum, download_filename, url, range = fail_lines[1].strip().split(',')
        assert acc == "GCA_000175535.1"
        assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
        assert moltype == "protein"
        assert download_filename == "GCA_000175535.1_protein.faa.gz"
        assert url == ""
        assert range == ""


def test_gbsketch_manifest(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    print(f"looking for path: {output}")

    idx = sourmash.load_file_as_index(output)
    manifest = sourmash_args.get_manifest(idx)
    assert len(manifest) == 3
    assert manifest._md5_set == set([ss1.md5sum(), ss2.md5sum(), ss3.md5sum()])
    for row in manifest.rows:
        print(row)
        if 'GCA_000175535.1' in row["name"]:
             assert row["md5"] == ss1.md5sum()
             assert row["n_hashes"] == 1047
        if "GCA_000961135.2" in row["name"]:
             if row["moltype"] == 'DNA':
                 assert row["md5"] == ss2.md5sum()
                 assert row["n_hashes"] == 1776
             else:
                 assert row["md5"] == ss3.md5sum()
                 assert row["n_hashes"] == 2596


def test_gbsketch_genomes_only(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--genomes-only',
                    '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 2
    for sig in sigs:
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            assert sig.md5sum() == ss2.md5sum()


def test_gbsketch_proteomes_only(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--proteomes-only',
                    '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        assert 'GCA_000961135.2' in sig.name
        assert sig.md5sum() == ss3.md5sum()


def test_gbsketch_genomes_only_via_params(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    captured = capfd.readouterr()

    assert len(sigs) == 2
    for sig in sigs:
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            assert sig.md5sum() == ss2.md5sum()
    assert 'No protein signature templates provided, and --keep-fasta is not set.' in captured.err
    assert 'Downloading and sketching genomes only.' in captured.err


def test_gbsketch_proteomes_only_via_params(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3',
                    '--checksum-fail', ch_fail,
                    '--param-str', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    print(runtmp.last_result.err)

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    captured = capfd.readouterr()

    assert len(sigs) == 1
    for sig in sigs:
        assert 'GCA_000961135.2' in sig.name
        assert sig.md5sum() == ss3.md5sum()
    assert 'No DNA signature templates provided, and --keep-fasta is not set.' in captured.err
    assert 'Downloading and sketching proteomes only.' in captured.err


def test_gbsketch_save_fastas(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--fastas', out_dir, '--keep-fasta',
                    '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.fna.gz', 'GCA_000961135.2_protein.faa.gz', 'GCA_000961135.2_genomic.fna.gz'])

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()


def test_gbsketch_save_fastas_no_overwrite(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    fa1 = runtmp.output('out_fastas/GCA_000961135.2_genomic.fna.gz')

    # first run to get one file downloaded 
    single_acc_csv = runtmp.output('single_acc.csv')
    with open(acc_csv, 'r') as inF, open(single_acc_csv, 'w') as outF:
        lines = inF.readlines()
        # only take the first acc for this test
        outF.write(lines[0])  # Header
        outF.write(lines[1])
    
    # Run the first gbsketch to download the fasta files
    runtmp.sourmash('scripts', 'gbsketch', single_acc_csv, '--download-only',
                    '--failed', failed, '-r', '3', '--fastas', out_dir, '--keep-fasta',
                    '--checksum-fail', ch_fail, '-g')
    assert os.path.exists(fa1)
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000961135.2_genomic.fna.gz'])

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '--download-only',
                    '--failed', failed, '-r', '3', '--fastas', out_dir, '--keep-fasta',
                    '--checksum-fail', ch_fail, '-g', '--no-overwrite-fasta')

    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    assert "Skipped 1 download(s) due to existing sketches and/or FASTA files." in captured.err

    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.fna.gz', 'GCA_000961135.2_genomic.fna.gz'])


def test_gbsketch_download_only(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '--download-only',
                    '--failed', failed, '-r', '3', '--fastas', out_dir, '--keep-fasta',
                    '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.fna.gz', 'GCA_000961135.2_protein.faa.gz', 'GCA_000961135.2_genomic.fna.gz'])
    captured = capfd.readouterr()
    print(captured)
    assert "Failed to send signatures: channel closed" not in captured.err


def test_gbsketch_bad_acc(runtmp):
    acc_csv = get_test_data('acc.csv')
    acc_mod = runtmp.output('acc_mod.csv')
    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        for line in lines:
            # if this acc exist in line, copy it and write an extra line with an invalid accession
            outF.write(line)
            print(line)
            if "GCA_000175535.1" in line:
                mod_line = line.replace('GCA_000175535.1', 'GCA_0001755559.1')  # add extra digit - should not be valid
                print(mod_line)
                outF.write(mod_line)

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '3', #'--fastas', output_fastas,
                    '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    # Open the failed file
    assert os.path.exists(failed)
    with open(failed, 'r') as acc_file:
        # Read the lines of the file
        lines = acc_file.readlines()
        # Check if the modified accession exists in the first column of any line
        for line in lines:
            print(line)
            if "GCA_0001755559.1" in line.split(',')[0]:
                assert True
                break
        else:
            assert False, "Modified accession not found"
    
    assert os.path.exists(failed)
    with open(ch_fail, 'r') as ch_file:
        # Read the lines of the file
        lines = ch_file.readlines()
        print("CHECKSUM FAILURES")
        print(lines)
        assert lines == ['accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason\n']

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()


def test_gbsketch_extra_column(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    acc_mod = runtmp.output('acc_mod.csv')

    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        for line in lines:
            outF.write(line.strip() + ',extra\n')

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    runtmp.sourmash('scripts', 'gbsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    print(f"looking for path: {output}")

    assert "WARNING: extra column 'extra' in CSV file. Ignoring." in captured.err


def test_gbsketch_bad_column_order(runtmp, capfd):
    acc_csv = get_test_data('acc-bad-header-order.csv')

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                        '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                        '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert not os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)

    assert "Invalid column names in CSV file." in captured.err


def test_gbsketch_missing_accfile(runtmp, capfd):
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No such file or directory" in captured.err


def test_gbsketch_empty_accfile(runtmp, capfd):
    acc_csv = get_test_data('acc1.csv')
    with open(acc_csv, 'w') as file:
        file.write('')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Invalid column names in CSV file." in captured.err


def test_gbsketch_bad_acc_fail(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    acc_mod = runtmp.output('acc_mod.csv')
    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        outF.write(lines[0])  # write the header line
        for line in lines:
            # if this acc exist in line, write bad acc instead
            if "GCA_000175535.1" in line:
                mod_line = line.replace('GCA_000175535.1', 'GCA_0001755559.1')  # add extra digit - should not be valid
                print(mod_line)
                outF.write(mod_line)
    
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000")
        
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Error: Failed to retrieve dehydrated download ZIP. Are your accessions valid? Exiting." in captured.err


def test_gbsketch_version_bug(runtmp):
    # test for bug where we didn't check version correctly
    acc_csv = get_test_data('acc-version.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000193795.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        assert sig.name == ss1.name == "GCA_000193795.2 Neisseria lactamica NS19 (b-proteobacteria) strain=NS19"
        assert sig.md5sum() == ss1.md5sum() == "7ead366dcfed8e8ab938f771459c4e94"


def test_gbsketch_cols_trailing_commas(runtmp, capfd):
    acc_csv = get_test_data('acc-cols.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error: CSV error: record 2 (line: 2, byte: 98): found record with 3 fields, but the previous record has 2 fields' in captured.err


def test_gbsketch_missing_output(runtmp):
    # no output sig zipfile provided but also not --download-only
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000")

    assert "Error: output signature zipfile is required if not using '--download-only'." in runtmp.last_result.err


def test_zip_file_permissions(runtmp):
    # Check permissions in the ZIP file
    import zipfile
    import stat
    
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out  # stdout should be empty

    with zipfile.ZipFile(output, 'r') as zip_ref:
        for zip_info in zip_ref.infolist():
            # The external_attr field contains the file permissions information.
            # By shifting right 16 bits (>> 16), we extract the file permissions.
            external_attr = zip_info.external_attr >> 16 
            permissions = stat.filemode(external_attr)
            print(f"File: {zip_info.filename}, Permissions: {permissions}")
            # check permissions are 644 (rw-r--r--)
            assert external_attr == 0o644


def test_gbsketch_protein_dayhoff_hp(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000961135.2.protein.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.dayhoff.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.hp.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=30, select_moltype='protein')
    ss2 = sourmash.load_one_signature(sig2, ksize=30, select_moltype='dayhoff')
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='hp')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3',
                    '--checksum-fail', ch_fail,
                    '--param-str',"protein,k=10,scaled=200",
                    '-p', "dayhoff,k=10,scaled=200",
                    '-p', "hp,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        assert sig.name == ss1.name
        if sig.minhash.moltype == 'protein':
            assert sig.md5sum() == ss1.md5sum()
        elif sig.minhash.moltype == 'dayhoff':
            assert sig.md5sum() == ss2.md5sum()
        elif sig.minhash.moltype == 'hp':
            assert sig.md5sum() == ss3.md5sum()
    assert os.path.exists(failed)
    with open(failed, 'r') as failF:
        fail_lines = failF.readlines()
        print(fail_lines)
        assert len(fail_lines) == 2
        assert fail_lines[0] == "accession,name,moltype,md5sum,download_filename,url,range\n"
        acc, name, moltype, md5sum, download_filename, url, range = fail_lines[1].strip().split(',')
        assert acc == "GCA_000175535.1"
        assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
        assert moltype == "protein"
        assert download_filename == "GCA_000175535.1_protein.faa.gz"
        assert url == ""
        assert range == ""


def test_gbsketch_simple_batched_single_acc(runtmp, capfd):
    # DNA and protein sigs for a single accession will be in separate zips now :shrug:
    acc_csv = get_test_data('acc.csv')
    acc1 = runtmp.output('acc1.csv')
    # open acc.csv with csv dictreader and keep accession= GCA_000961135.2 line
    with open(acc_csv, 'r') as inF, open(acc1, 'w', newline='') as outF:
        r = csv.DictReader(inF)
        w = csv.DictWriter(outF, fieldnames=r.fieldnames)
        w.writeheader()
        for row in r:
            if row['accession'] == "GCA_000961135.2":
                w.writerow(row)

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')

    sig1 = get_test_data('GCA_000961135.2.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc1, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)

    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss1.name, ss2.md5sum(), ss2.minhash.moltype), # ss1 name b/c of how it's written in acc.csv
    }

    # Collect the actual signature information from all the output files
    all_siginfo = set()
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    for sig in sigs:
        all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))
    idx2 = sourmash.load_file_as_index(out2)
    sigs2 = list(idx2.signatures())
    for sig in sigs2:
        all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_batched_multiple(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '2')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)

    assert "finished batch 1" in captured.err
    assert "finished batch 2" in captured.err
    print(captured.out)
    print(runtmp.last_result.err)
    batch_base = output.split('.zip')[0]
    print(batch_base)
    assert f"Sigs in '{batch_base}.1.zip', etc" in runtmp.last_result.err

    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss2.name, ss2.md5sum(), ss2.minhash.moltype),
        (ss2.name, ss3.md5sum(), ss3.minhash.moltype), # ss2 name b/c of how it's written in acc.csv
    }

    # Collect the actual signature information from all the output files
    all_siginfo = set()
    for out_file in [out1, out2]:
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_batch_restart(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')
    out3 = runtmp.output('simple.3.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)
    # why does this need ksize =30 and not ksize = 10!???
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)

    assert "Skipped 1 download(s) due to missing download URLs." in captured.err

    # # we created this one with sig cat
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    assert len(sigs) == 2
    for sig in sigs:
        assert sig.name == ss2.name
        assert sig.md5sum() in [ss2.md5sum(), ss3.md5sum()]

    # # these were created with gbsketch
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
    }

    # Collect actual signature information from gbsketch zip batches
    all_siginfo = set()
    for out_file in [out2, out3]:
        print( f"Loading signatures from {out_file}...")
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_batch_restart_sig_zip(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.sig.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.sig.zip')
    out2 = runtmp.output('simple.2.sig.zip')
    out3 = runtmp.output('simple.3.sig.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)
    # why does this need ksize =30 and not ksize = 10!???
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)

    # # we created this one with sig cat
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    assert len(sigs) == 2
    for sig in sigs:
        assert sig.name == ss2.name
        assert sig.md5sum() in [ss2.md5sum(), ss3.md5sum()]

    # # these were created with gbsketch
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
    }

    # Collect actual signature information from gbsketch zip batches
    all_siginfo = set()
    for out_file in [out2, out3]:
        print( f"Loading signatures from {out_file}...")
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_batch_restart_incomplete_sig_zip(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.sig.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.sig.zip')
    out2 = runtmp.output('simple.2.sig.zip')
    out2_incomplete = runtmp.output('simple.2.sig.zip.incomplete')
    out3 = runtmp.output('simple.3.sig.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)
    # why does this need ksize =30 and not ksize = 10!???
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)
    with open(out2_incomplete, 'w') as f:
        pass  # create an empty out2 to simulate an incomplete batch

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)
    print(captured.out)

    assert "Removing incomplete zip file: 'simple.2.sig.zip.incomplete'" in captured.err

    # # we created this one with sig cat
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    assert len(sigs) == 2
    for sig in sigs:
        assert sig.name == ss2.name
        assert sig.md5sum() in [ss2.md5sum(), ss3.md5sum()]

    # # these were created with gbsketch
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
    }

    # Collect actual signature information from gbsketch zip batches
    all_siginfo = set()
    for out_file in [out2, out3]:
        print( f"Loading signatures from {out_file}...")
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_batch_restart_skipcount(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')
    out3 = runtmp.output('simple.3.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)

    assert "Found 1 existing valid zip batch(es). Starting new sig writing at batch 2" in captured.err
    assert "Skipped 1 download(s) due to existing sketches and/or FASTA files." in captured.err

    # # we created this one with sig cat
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    assert len(sigs) == 2
    for sig in sigs:
        assert sig.name == ss2.name
        assert sig.md5sum() in [ss2.md5sum(), ss3.md5sum()]

    # # these were created with gbsketch
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
    }

    # Collect actual signature information from gbsketch zip batches
    all_siginfo = set()
    for out_file in [out2]:
        print( f"Loading signatures from {out_file}...")
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_negative_batch_size(runtmp):
    # negative int provided for batch size
    acc_csv = runtmp.output('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv,
                    '--failed', failed, '-r', '3', '--batch-size', '-2',
                    '--param-str', "dna,k=31,scaled=1000")

    assert "Batch size cannot be negative (input value: -2)" in runtmp.last_result.err


def test_gbsketch_simple_batch_restart_with_incomplete_zip(runtmp, capfd):
    # test restart with complete + incomplete zipfile batches
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')
    out3 = runtmp.output('simple.3.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)

    # Create an invalid zip file for out2
    with open(out2, 'wb') as f:
        f.write(b"This is not a valid zip file!")

    assert os.path.exists(out2)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)  # Should be overwritten
    assert os.path.exists(out3)
    assert not os.path.exists(output)  # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)
    assert f"Warning: Failed to load zip file '{out2}'; skipping." in captured.err

    # we created this one with sig cat
    idx = sourmash.load_file_as_index(out1)
    sigs = list(idx.signatures())
    assert len(sigs) == 2
    for sig in sigs:
        assert sig.name == ss2.name
        assert ss2.md5sum() in [ss2.md5sum(), ss3.md5sum()]

    # these were created with gbsketch (out2 should have been overwritten)
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
    }

    # Collect actual signature information from gbsketch zip batches
    all_siginfo = set()
    for out_file in [out2, out3]:
        # this would fail if out2 were not overwritten with a valid sig zip
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Assert that all expected signatures are found (ignoring order)
    assert all_siginfo == expected_siginfo


def test_gbsketch_bad_param_str(runtmp, capfd):
    # negative int provided for batch size
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('ch.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, "-o", output,
                    '--failed', failed, '-r', '3', "--checksum-fail", ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,protein")
        
    captured = capfd.readouterr()
    print(captured)

    assert "Failed to parse params string: Conflicting moltype settings in param string: 'dna' and 'protein'" in captured.err


def test_gbsketch_overwrite(runtmp, capfd):
    # test restart with complete + incomplete zipfile batches
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')
    out1 = runtmp.output('simple.1.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # ss3 = sourmash.load_one_signature(sig2, ksize=21)
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # run the workflow once - write all to single output
    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200")
    assert os.path.exists(output)
    captured = capfd.readouterr()
    print(captured.err)
    expected_siginfo = {
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
        (ss2.name, ss2.md5sum(), ss2.minhash.moltype),
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
    }

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    all_siginfo = set()
    for sig in sigs:
        all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    assert all_siginfo == expected_siginfo

    # now, try running again - providing same output file.
    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200")

    # check that sigs can still be read
    assert os.path.exists(output)
    assert not os.path.exists(out1)
    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    all_siginfo = set()
    for sig in sigs:
        all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))
    assert all_siginfo == expected_siginfo


def test_gbsketch_simple_skipmer(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '3', '--checksum-fail', ch_fail,
                    '--param-str', "skipm2n3,k=21,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.err)
    print(f"looking for path: {output}")

    # read the file with python and check sigs
    import zipfile, gzip, json

    with zipfile.ZipFile(output, "r") as zf:
        # Check the manifest exists
        assert "SOURMASH-MANIFEST.csv" in zf.namelist()

        expected_signatures = [
            {
                "name": "GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454",
                "ksize": 21,
                "scaled": 1000,
                "moltype": "skipm2n3",
                "md5sum": "5745400ada0c3a27ddf1e8d5d1a46b7a",
            },
            {
                "name": "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14",
                "ksize": 21,
                "scaled": 1000,
                "moltype": "skipm2n3",
                "md5sum": "2e37ca0bb9228bc3f5e1337e8535ab26",
            },
        ]
        expected_signatures_dict = {exp["md5sum"]: exp for exp in expected_signatures}

        # Read and parse the manifest
        with zf.open("SOURMASH-MANIFEST.csv") as manifest_file:
            manifest_data = manifest_file.read().decode("utf-8").splitlines()
            manifest_data = [line for line in manifest_data if not line.startswith("#")]
            manifest_reader = csv.DictReader(manifest_data)

            for row in manifest_reader:
                if row["moltype"] == "skipm2n3":
                    print("Manifest Row:", row)

                    # Validate row fields
                    md5sum = row["md5"]
                    assert (
                        md5sum in expected_signatures_dict
                    ), f"Unexpected md5sum: {md5sum}"
                    expected = expected_signatures_dict[md5sum]
                    assert (
                        row["name"] == expected["name"]
                    ), f"Name mismatch: {row['name']}"
                    assert (
                        int(row["ksize"]) == expected["ksize"]
                    ), f"Ksize mismatch: {row['ksize']}"
                    assert (
                        row["moltype"] == expected["moltype"]
                    ), f"Moltype mismatch: {row['moltype']}"

                    sig_path = row["internal_location"]
                    assert sig_path.startswith("signatures/")

                    # Extract and read the signature file
                    with zf.open(sig_path) as sig_gz:
                        with gzip.open(sig_gz, "rt") as sig_file:
                            sig_contents = json.load(sig_file)
                            print("Signature Contents:", sig_contents)

                            # Validate signature contents
                            sig_data = sig_contents[0]
                            print(sig_data)
                            siginfo = sig_data["signatures"][0]
                            assert (
                                siginfo["md5sum"] == md5sum
                            ), f"MD5 mismatch: {siginfo['md5sum']}"
                            assert (
                                siginfo["ksize"] == expected["ksize"]
                            ), f"Ksize mismatch: {siginfo['ksize']}"
                            assert (
                                siginfo["molecule"] == expected["moltype"]
                            ), f"Moltype mismatch: {siginfo['molecule']}"


def test_gbsketch_max_n_downloads(runtmp, capfd):
    #check that we can use 30 simultaneous downloads
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    assert  os.environ["NCBI_API_KEY"] == ""

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-n', '30', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)

    assert "using 30 simultaneous downloads, 3 retries" in captured.out
    assert "Successfully downloaded and parsed dehydrated zipfile. Now processing accessions." in captured.err


def test_gbsketch_verbose(runtmp, capfd):
    #test verbose reporting
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    assert  os.environ["NCBI_API_KEY"] == ""

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-n', '30', '--checksum-fail', ch_fail, "--verbose",
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)

    assert "Starting download 1/3 (33%) - accession: 'GCA_000961135.2', moltype: DNA" in captured.out
    assert "Starting download 2/3 (67%) - accession: 'GCA_000961135.2', moltype: protein" in captured.out
    assert "Starting download 3/3 (100%) - accession: 'GCA_000175535.1', moltype: DNA" in captured.out


def test_gbsketch_from_gbsketch_failed(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(failed)
    with open(failed, 'r') as failF:
        fail_lines = failF.readlines()
        assert len(fail_lines) == 2
        assert fail_lines[0] == "accession,name,moltype,md5sum,download_filename,url,range\n"
        acc, name, moltype, md5sum, download_filename, url, range = fail_lines[1].strip().split(',')
        assert acc == "GCA_000175535.1"
        assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
        assert moltype == "protein"
        assert download_filename == "GCA_000175535.1_protein.faa.gz"
        assert url == ""
        assert range == ""
    assert not runtmp.last_result.out # stdout should be empty

    out2 = runtmp.output('failed-retry.zip')
    fail2 = runtmp.output('fail2.csv')

    # since the protein file doesn't exist at NCBI, we won't be able to find any links in the downloaded dehydrated zip.
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', failed, '-o', out2,
                    '--failed', fail2, '-r', '1',
                    '-p', "protein,k=10,scaled=200")
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Are your accessions valid?" in captured.err


def test_gbsketch_write_urlsketch_csv(runtmp, capfd):
    #test verbose reporting
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')
    urlsketch_csv = runtmp.output(f'{acc_csv}.urlsketch.csv')

    assert  os.environ["NCBI_API_KEY"] == ""

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-n', '30', '--checksum-fail', ch_fail, "--verbose",
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200",
                    '--write-urlsketch-csv')
    print(runtmp.last_result.out)
    print(runtmp.last_result.err)
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)

    assert os.path.exists(urlsketch_csv), f"URL sketch CSV not found: {urlsketch_csv}"
    with open(urlsketch_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        assert len(rows) == 3  # 3 entries (excluding the header)

        # Check one full line
        first_row = rows[0]
        assert first_row["accession"] == "GCA_000961135.2"
        assert first_row["name"] == "GCA_000961135.2 Candidatus Aramenus sulfurataquae isolate AZ1-454"
        assert first_row["moltype"] == "DNA"
        assert first_row["download_filename"] == "GCA_000961135.2_genomic.fna.gz"
        assert first_row["url"] == "https://api.ncbi.nlm.nih.gov/datasets/fetch_h/R2V0UmVtb3RlRGF0YWZpbGU/eNqTyuXKzktOytROKymw0tdPT83Lz00t1k_MydF3d3bUNzAw0Lc0M9Q3NDYF8eOBfCAXyNMzincM9gWzy4zwSMWDTcxM1kvLS9RLrzJgtGAEADqOH0U"

        second_row = rows[1]
        assert second_row["accession"] == "GCA_000961135.2"
        assert second_row["moltype"] == "protein"
        assert second_row["url"].startswith('https://api.ncbi.nlm.nih.gov/datasets/fetch')
        third_row = rows[2]
        assert third_row["accession"] == "GCA_000175535.1"
        assert third_row['moltype'] == "DNA"
        assert third_row["url"].startswith('https://api.ncbi.nlm.nih.gov/datasets/fetch')
        # Print all rows for debugging
        for row in rows:
            print(row)

    assert "Starting download 1/3 (33%) - accession: 'GCA_000961135.2', moltype: DNA" in captured.out
    assert "Starting download 2/3 (67%) - accession: 'GCA_000961135.2', moltype: protein" in captured.out
    assert "Starting download 3/3 (100%) - accession: 'GCA_000175535.1', moltype: DNA" in captured.out
