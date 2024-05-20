"""
Tests for urlsketch
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
        runtmp.sourmash('scripts', 'urlsketch')

    assert 'usage:  urlsketch' in runtmp.last_result.err


def test_urlsketch_simple(runtmp):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

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
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url = line.strip().split(',')
            assert acc == "GCA_000175535.1"
            assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
            assert moltype == "protein"
            assert download_filename == "GCA_000175535.1_protein.faa.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_protein.faa.gz"

def test_urlsketch_save_fastas(runtmp):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')


    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fasta',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.urlsketch.fna.gz', 'GCA_000961135.2_protein.urlsketch.faa.gz', 'GCA_000961135.2_genomic.urlsketch.fna.gz'])

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


def test_urlsketch_download_only(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')


    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '--download-only',
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fasta',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.urlsketch.fna.gz', 'GCA_000961135.2_protein.urlsketch.faa.gz', 'GCA_000961135.2_genomic.urlsketch.fna.gz'])
    captured = capfd.readouterr()
    assert "Failed to send signatures: channel closed" not in captured.err