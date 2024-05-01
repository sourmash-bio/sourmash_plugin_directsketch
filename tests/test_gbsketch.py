"""
Tests for gbsketch
"""
import os
import pytest

import sourmash
import sourmash_tst_utils as utils
from sourmash_tst_utils import SourmashCommandFailed



def get_test_data(filename):
    thisdir = os.path.dirname(__file__)
    return os.path.join(thisdir, 'test-data', filename)

def test_installed(runtmp):
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch')

    assert 'usage:  gbsketch' in runtmp.last_result.err


def test_gbsketch_simple(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175555.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175555.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()

def test_gbsketch_genomes_only(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175555.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--genomes-only',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 2
    for sig in sigs:
        if 'GCA_000175555.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            assert sig.md5sum() == ss2.md5sum()


def test_gbsketch_proteomes_only(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--proteomes-only',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        assert 'GCA_000961135.2' in sig.name
        assert sig.md5sum() == ss3.md5sum()


def test_gbsketch_save_fastas(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')


    sig1 = get_test_data('GCA_000175555.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fastas',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175555.1_genomic.fna.gz', 'GCA_000961135.2_protein.faa.gz', 'GCA_000961135.2_genomic.fna.gz'])

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175555.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()

def test_gbsketch_download_only(runtmp):
    acc_csv = get_test_data('acc.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')


    sig1 = get_test_data('GCA_000175555.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output, '--download-only',
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fastas',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output) # would be better if this didn't exist
    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175555.1_genomic.fna.gz', 'GCA_000961135.2_protein.faa.gz', 'GCA_000961135.2_genomic.fna.gz'])
   

def test_gbsketch_bad_acc(runtmp):
    acc_csv = get_test_data('acc.csv')
    acc_mod = runtmp.output('acc_mod.csv')
    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        for line in lines:
            # if this acc exist in line, copy it and write an extra line with an invalid accession
            outF.write(line)
            print(line)
            if "GCA_000175555.1" in line:
                mod_line = line.replace('GCA_000175555.1', 'GCA_0001755559.1')  # add extra digit - should not be valid
                print(mod_line)
                outF.write(mod_line)

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175555.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'gbsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1', #'--fastas', output_fastas,
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

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 3
    for sig in sigs:
        if 'GCA_000175555.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()
    

def test_gbsketch_missing_accfile(runtmp, capfd):
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
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

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No accessions to download and sketch." in captured.err


def test_gbsketch_bad_acc_fail(runtmp, capfd):
    acc_csv = get_test_data('acc.csv')
    acc_mod = runtmp.output('acc_mod.csv')
    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        outF.write(lines[0])  # write the header line
        for line in lines:
            # if this acc exist in line, copy it and write
            if "GCA_000175555.1" in line:
                mod_line = line.replace('GCA_000175555.1', 'GCA_0001755559.1')  # add extra digit - should not be valid
                print(mod_line)
                outF.write(mod_line)
    
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'gbsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")
        
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Error: No signatures written." in captured.err