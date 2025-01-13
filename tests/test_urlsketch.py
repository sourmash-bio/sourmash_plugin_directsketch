"""
Tests for urlsketch
"""
import os
import pytest
import gzip
import screed

import csv
import sourmash
from sourmash import sourmash_args
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
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "GCA_000175535.1"
            assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
            assert moltype == "protein"
            assert download_filename == "GCA_000175535.1_protein.faa.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_protein.faa.gz"
            assert range == ""


def test_urlsketch_manifest(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
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


def test_urlsketch_save_fastas_no_append_across_runs(runtmp):
    # make sure we overwrite files on subsequent runs (not append to existing)
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')

    # run once
    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fasta',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    # check out fastas exist
    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    fa_files = os.listdir(out_dir)
    assert set(fa_files) == set(['GCA_000175535.1_genomic.urlsketch.fna.gz', 'GCA_000961135.2_protein.urlsketch.faa.gz', 'GCA_000961135.2_genomic.urlsketch.fna.gz'])

    # Get the file size for each file
    fsizes = set()
    for fa_file in fa_files:
        file_path = os.path.join(out_dir, fa_file)
        file_size = os.path.getsize(file_path)
        print(f"File: {fa_file}, Size: {file_size} bytes")
        fsizes.add(file_size)

    # run a second time
    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--fastas', out_dir, '--keep-fasta',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")

    fa_files2 = os.listdir(out_dir)
    assert set(fa_files2) == set(['GCA_000175535.1_genomic.urlsketch.fna.gz', 'GCA_000961135.2_protein.urlsketch.faa.gz', 'GCA_000961135.2_genomic.urlsketch.fna.gz'])
    for fa_file in fa_files2:
        file_path = os.path.join(out_dir, fa_file)
        file_size = os.path.getsize(file_path)
        print(f"File: {fa_file}, Size: {file_size} bytes")
        assert file_size in fsizes


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


def test_urlsketch_bad_acc(runtmp):
    acc_csv = get_test_data('acc-url.csv')
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

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    # why does this need ksize =30 and not ksize = 10!???
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
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
        if 'GCA_000175535.1' in sig.name:
            assert sig.name == ss1.name
            assert sig.md5sum() == ss1.md5sum()
        elif 'GCA_000961135.2' in sig.name:
            assert sig.name == ss2.name
            if sig.minhash.moltype == 'DNA':
                assert sig.md5sum() == ss2.md5sum()
            else:
                assert sig.md5sum() == ss3.md5sum()


def test_urlsketch_missing_accfile(runtmp, capfd):
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: No such file or directory" in captured.err


def test_urlsketch_empty_accfile(runtmp, capfd):
    acc_csv = get_test_data('acc1.csv')
    with open(acc_csv, 'w') as file:
        file.write('')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200")
        
    captured = capfd.readouterr()
    print(captured.err)
    assert 'Error: Invalid column names in CSV file. Columns should be: ["accession", "name", "moltype", "md5sum", "download_filename", "url", "range"]' in captured.err


def test_urlsketch_bad_acc_fail(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    acc_mod = runtmp.output('acc_mod.csv')
    with open(acc_csv, 'r') as inF, open(acc_mod, 'w') as outF:
        lines = inF.readlines()
        outF.write(lines[0])  # write the header line
        for line in lines:
            # if this acc exist in line, copy it and write
            if "GCA_000175535.1" in line:
                mod_line = line.replace('GCA_000175535.1', 'GCA_0001755559.1')  # add extra digit - should not be valid
                print(mod_line)
                outF.write(mod_line)
    
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")
        
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Error: No signatures written, exiting." in captured.err


def test_urlsketch_missing_output(runtmp):
    # no output sig zipfile provided but also not --download-only
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_csv,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")

    assert "Error: output signature zipfile is required if not using '--download-only'." in runtmp.last_result.err


def test_urlsketch_from_gbsketch_failed(runtmp, capfd):
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
        assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_protein.faa.gz"
        assert range == ""
    assert not runtmp.last_result.out # stdout should be empty

    out2 = runtmp.output('failed-retry.zip')
    fail2 = runtmp.output('fail2.csv')

    with pytest.raises(utils.SourmashCommandFailed):

        runtmp.sourmash('scripts', 'urlsketch', failed, '-o', out2,
                    '--failed', fail2, '-r', '1',
                    '-p', "protein,k=10,scaled=200")
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)
    assert "Error: No signatures written, exiting." in captured.err

    # since no protein file exists, fail2 should just be the same as failed
    assert os.path.exists(fail2)
    with open(fail2, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "GCA_000175535.1"
            assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
            assert moltype == "protein"
            assert download_filename == "GCA_000175535.1_protein.faa.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_protein.faa.gz"
            assert range == ""


def test_zip_file_permissions(runtmp):
    # Check permissions in the ZIP file
    import zipfile
    import stat
    
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
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


def test_urlsketch_protein_dayhoff_hp(runtmp):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000961135.2.protein.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.dayhoff.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.hp.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=30, select_moltype='protein')
    ss2 = sourmash.load_one_signature(sig2, ksize=30, select_moltype='dayhoff')
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='hp')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
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
        assert len(fail_lines) == 1
        assert fail_lines[0] == "accession,name,moltype,md5sum,download_filename,url,range\n"


def test_urlsketch_md5sum_mismatch_checksum_file(runtmp, capfd):
    acc_csv = get_test_data('acc-url-md5sum.csv')

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    sig1 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        assert sig.name == ss1.name
        assert sig.md5sum() == ss1.md5sum()

    assert os.path.exists(ch_fail)
    with open(ch_fail, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum_url, download_filename, url, expected_md5, reason= line.strip().split(',')
            assert acc == "GCA_000175535.1"
            assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
            assert moltype == "DNA"
            assert md5sum_url == ""
            assert expected_md5 == "b1234567"
            assert download_filename == "GCA_000175535.1_genomic.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert reason == "MD5 hash does not match. Expected: 'b1234567'; Found: 'a1a8f1c6dc56999c73fe298871c963d1'"


def test_urlsketch_md5sum_mismatch_no_checksum_file(runtmp, capfd):
    acc_csv = get_test_data('acc-url-md5sum.csv')

    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000961135.2.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    captured = capfd.readouterr()
    print(captured.out)
    print(captured.err)

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        assert sig.name == ss1.name
        assert sig.md5sum() == ss1.md5sum()

    assert os.path.exists(failed)
    with open(failed, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range= line.strip().split(',')
            assert acc == "GCA_000175535.1"
            assert name == "GCA_000175535.1 Chlamydia muridarum MopnTet14 (agent of mouse pneumonitis) strain=MopnTet14"
            assert moltype == "DNA"
            assert md5sum == "b1234567"
            assert download_filename == "GCA_000175535.1_genomic.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == ""


def test_urlsketch_simple_batched(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')
    
    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')
    out3 = runtmp.output('simple.3.zip')
    out4 = runtmp.output('simple.4.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(out4)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)
    assert "finished batch 1" in captured.err
    assert "finished batch 2" in captured.err
    assert "finished batch 3" in captured.err
    print(captured.out)
    print(runtmp.last_result.err)
    batch_base = output.split('.zip')[0]
    print(batch_base)
    assert f"Sigs in '{batch_base}.1.zip', etc" in runtmp.last_result.err

    expected_siginfo = {
            (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
            (ss2.name, ss2.md5sum(), ss2.minhash.moltype),
            (ss3.name, ss3.md5sum(), ss3.minhash.moltype)
        }
    # Collect all signatures from the output zip files
    all_sigs = []

    for out_file in [out1, out2, out3]:
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        assert len(sigs) == 1  # We expect exactly 1 signature per batch
        all_sigs.append(sigs[0])
    
    loaded_signatures = {(sig.name, sig.md5sum(), sig.minhash.moltype) for sig in all_sigs}
    assert loaded_signatures == expected_siginfo, f"Loaded sigs: {loaded_signatures}, expected: {expected_siginfo}"


def test_urlsketch_simple_batch_restart(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('simple.1.zip')
    out2 = runtmp.output('simple.2.zip')
    out3 = runtmp.output('simple.3.zip')
    out4 = runtmp.output('simple.4.zip')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    sig3 = get_test_data('GCA_000961135.2.protein.sig.gz')
    ss1 = sourmash.load_one_signature(sig1, ksize=31)
    ss2 = sourmash.load_one_signature(sig2, ksize=31)
    ss3 = sourmash.load_one_signature(sig2, ksize=21)
    ss4 = sourmash.load_one_signature(sig3, ksize=30, select_moltype='protein')

    # first, cat sig2 into an output file that will trick gbsketch into thinking it's a prior batch
    # need to actually rename it first, so it will match sig that would have been written

    runtmp.sourmash('sig', 'cat', sig2, '-o', out1)
    assert os.path.exists(out1)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '5', '-n', "1", '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(out4)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)
  
    expected_siginfo = {
        (ss2.name, ss2.md5sum(), ss2.minhash.moltype),
        (ss2.name, ss3.md5sum(), ss3.minhash.moltype), # ss2 name b/c thats how it is in acc-url.csv
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
    }

    all_siginfo = set()
    for out_file in [out1, out2, out3]:
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Verify that the loaded signatures match the expected signatures, order-independent
    assert all_siginfo == expected_siginfo, f"Loaded sigs: {all_siginfo}, expected: {expected_siginfo}"


def test_urlsketch_negative_batch_size(runtmp):
    # negative int provided for batch size
    acc_csv = runtmp.output('acc1.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_csv,
                    '--failed', failed, '-r', '1', '--batch-size', '-2',
                    '--param-str', "dna,k=31,scaled=1000")

    assert "Batch size cannot be negative (input value: -2)" in runtmp.last_result.err


def test_urlsketch_simple_batch_restart_with_incomplete_zip(runtmp, capfd):
    # test restart with complete + incomplete zipfile batches
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('restart.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    out1 = runtmp.output('restart.1.zip')
    out2 = runtmp.output('restart.2.zip')
    out3 = runtmp.output('restart.3.zip')
    out4 = runtmp.output('restart.4.zip')

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

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '5', '-n', "1", '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=1000,abund", '-p', "protein,k=10,scaled=200",
                    '--batch-size', '1')

    assert os.path.exists(out1)
    assert os.path.exists(out2)
    assert os.path.exists(out3)
    assert not os.path.exists(out4)
    assert not os.path.exists(output) # for now, orig output file should be empty.
    captured = capfd.readouterr()
    print(captured.err)
    assert f"Warning: Failed to load zip file '{out2}'; skipping." in captured.err

    expected_siginfo = {
        (ss2.name, ss2.md5sum(), ss2.minhash.moltype),
        (ss2.name, ss3.md5sum(), ss3.minhash.moltype), # ss2 name b/c thats how it is in acc-url.csv
        (ss4.name, ss4.md5sum(), ss4.minhash.moltype),
        (ss1.name, ss1.md5sum(), ss1.minhash.moltype),
    }

    all_siginfo = set()
    for out_file in [out1, out2, out3]:
        idx = sourmash.load_file_as_index(out_file)
        sigs = list(idx.signatures())
        for sig in sigs:
            all_siginfo.add((sig.name, sig.md5sum(), sig.minhash.moltype))

    # Verify that the loaded signatures match the expected signatures, order-independent
    assert all_siginfo == expected_siginfo, f"Loaded sigs: {all_siginfo}, expected: {expected_siginfo}"


def test_urlsketch_simple_skipmer(runtmp, capfd):
    acc_csv = get_test_data('acc-url.csv')
    output = runtmp.output('simple.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('checksum_dl_failed.csv')

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--checksum-fail', ch_fail,
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


def test_urlsketch_simple_merged(runtmp):
    acc_csv = get_test_data('acc-merged.csv')
    output = runtmp.output('merged.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    merged_sig = runtmp.output("sigmerge.zip")

    # create merged signature
    runtmp.sourmash("sig", "merge", "-k", "31", sig1, sig2, "--set-name", "both name", '-o', merged_sig)
    msigidx = sourmash.load_file_as_index(merged_sig)
    msig = list(msigidx.signatures())[0]
    print(msig.name)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    sig = sigs[0]
    assert sig.name == msig.name == "both name"
    print(msig.md5sum())
    assert sig.md5sum() == msig.md5sum()
    assert sig.minhash.moltype == msig.minhash.moltype == "DNA"
    assert os.path.exists(failed)


def test_urlsketch_simple_merged_with_md5sums(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    output = runtmp.output('merged.zip')
    failed = runtmp.output('failed.csv')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    merged_sig = runtmp.output("sigmerge.zip")

    # create merged signature
    runtmp.sourmash("sig", "merge", "-k", "31", sig1, sig2, "--set-name", "both name", '-o', merged_sig)
    msigidx = sourmash.load_file_as_index(merged_sig)
    msig = list(msigidx.signatures())[0]
    print(msig.name)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    sig = sigs[0]
    assert sig.name == msig.name == "both name"
    print(msig.md5sum())
    assert sig.md5sum() == msig.md5sum()
    assert sig.minhash.moltype == msig.minhash.moltype == "DNA"
    assert os.path.exists(failed)


def test_urlsketch_simple_merged_keep_fasta(runtmp):
    acc_csv = get_test_data('acc-merged.csv')
    output = runtmp.output('merged.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    merged_sig = runtmp.output("sigmerge.zip")

    # create merged signature
    runtmp.sourmash("sig", "merge", "-k", "31", sig1, sig2, "--set-name", "both name", '-o', merged_sig)
    msigidx = sourmash.load_file_as_index(merged_sig)
    msig = list(msigidx.signatures())[0]
    print(msig.name)

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--keep-fasta',
                    '--fastas', out_dir,
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    # check fasta files are present
    fa_files = os.listdir(out_dir)
    print(fa_files)
    assert fa_files == ['both.urlsketch.fna.gz']

    # check fasta files have records from both entries
    n_expected_records = 104
    n_records = 0
    # check one record from each
    expected_names = ["ACUJ01000001.1 Chlamydia muridarum MopnTet14 chromosome, whole genome shotgun sequence",
                      "JZWS02000016.1 MAG: Candidatus Aramenus sulfurataquae isolate AZ1-454 NODE_87_length_15535_cov_30.701232, whole genome shotgun sequence"]
    rec_names = []
    with screed.open(os.path.join(out_dir, fa_files[0])) as inF:
        for rec in inF:
            n_records +=1
            rec_names.append(rec.name)

    assert n_records == n_expected_records
    assert all(n in rec_names for n in expected_names)

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    sig = sigs[0]
    assert sig.name == msig.name == "both name"
    print(msig.md5sum())
    assert sig.md5sum() == msig.md5sum()
    assert sig.minhash.moltype == msig.minhash.moltype == "DNA"
    assert os.path.exists(failed)


def test_urlsketch_simple_merged_keep_fasta_path_in_filename(runtmp):
    acc_csv = get_test_data('acc-merged.csv')
    mod_csv = runtmp.output('acc-merged-filepath.csv')
    output = runtmp.output('merged.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')

    # open acc-merged.csv and prepend "/unavailable-path/subdir/" to the "download_filename" column
    with open(acc_csv, 'r') as infile, open(mod_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row['download_filename'] = f"unavailable-path/subdir/{row['download_filename']}"
            writer.writerow(row)

    sig1 = get_test_data('GCA_000175535.1.sig.gz')
    sig2 = get_test_data('GCA_000961135.2.sig.gz')
    merged_sig = runtmp.output("sigmerge.zip")

    # create merged signature
    runtmp.sourmash("sig", "merge", "-k", "31", sig1, sig2, "--set-name", "both name", '-o', merged_sig)
    msigidx = sourmash.load_file_as_index(merged_sig)
    msig = list(msigidx.signatures())[0]
    print(msig.name)

    runtmp.sourmash('scripts', 'urlsketch', mod_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--keep-fasta',
                    '--fastas', out_dir,
                    '--param-str', "dna,k=31,scaled=1000")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    # check fasta files are present
    fa_files = []
    for root, dirs, files in os.walk(out_dir):
        for file in files:
            if file.endswith('fna.gz'):
                fa_files.append(os.path.relpath(os.path.join(root, file), out_dir))
    print(fa_files)
    assert fa_files == ['unavailable-path/subdir/both.urlsketch.fna.gz']


def test_urlsketch_simple_merged_incorrect_md5sum_checksum_failure(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    mod_csv = runtmp.output('acc-merged_incorrect_md5.csv')
    output = runtmp.output('merged.zip')
    failed = runtmp.output('failed.csv')
    ch_failed = runtmp.output('ch-failed.csv')
    out_dir = runtmp.output('out_fastas')

    # open file and write incorrect md5sum
    with open(acc_csv, 'r') as infile, open(mod_csv, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row['md5sum'] = row['md5sum'][2:] # take off first digit from first md5sum
            print(row)
            writer.writerow(row)

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', mod_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--keep-fasta',
                    '--fastas', out_dir, '--checksum-fail', ch_failed,
                    '--param-str', "dna,k=31,scaled=1000")

    assert not os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    # check failure file
    assert os.path.exists(ch_failed)
    with open(ch_failed, 'r') as failF:
        header = next(failF).strip()
        print(header)
        assert header == "accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum_url, download_filename, url, expected_md5sum, reason = line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert download_filename == "both.urlsketch.fna.gz"
            assert expected_md5sum == "b9fb20c51f0552b87db5d44d5d4566"
            assert reason == "MD5 hash does not match. Expected: 'b9fb20c51f0552b87db5d44d5d4566'; Found: '47b9fb20c51f0552b87db5d44d5d4566'"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_genomic.fna.gz"
    with open(failed, 'r') as fails:
        header = next(fails).strip()
        print(header)
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in fails:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert md5sum == "b9fb20c51f0552b87db5d44d5d4566;a1a8f1c6dc56999c73fe298871c963d1"
            assert download_filename == "both.urlsketch.fna.gz"
            assert url ==  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == ""


def test_urlsketch_with_range(runtmp):
    acc_csv = get_test_data('acc-url-range.csv')
    subseqs = get_test_data('subseqs.zip')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')

    # open subseq sigs
    idx = sourmash.load_file_as_index(subseqs)
    siglist = list(idx.signatures())
    ss1 = siglist[0]
    ss2 = siglist[1]

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 2
    for sig in sigs:
        ident = sig.name.split(' ')[0]
        assert ident in ["GCA_000175535.1_first50kb", "GCA_000175535.1_second50kb"]
        print(ident)
        if ident == "GCA_000175535.1_first50kb":
            assert sig.md5sum() == ss1.md5sum()
        if ident == "GCA_000175535.1_second50kb":
            assert sig.md5sum() == ss2.md5sum()
    assert os.path.exists(failed)


def test_urlsketch_with_range_keep_fasta(runtmp):
    acc_csv = get_test_data('acc-url-range.csv')
    subseqs = get_test_data('subseqs.zip')
    first50kb = get_test_data('GCA_000175535.1_ASM17553v1_genomic.1-50000.fna.gz')
    second50kb = get_test_data('GCA_000175535.1_ASM17553v1_genomic.50000-100000.fna.gz')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')
    out_dir = runtmp.output('out_fastas')

    # open subseq sigs
    idx = sourmash.load_file_as_index(subseqs)
    siglist = list(idx.signatures())
    ss1 = siglist[0]
    ss2 = siglist[1]

    runtmp.sourmash('scripts', 'urlsketch', acc_csv, '-o', output,
                    '--failed', failed, '-r', '1', '--keep-fasta',
                    '--fastas', out_dir,
                    '--param-str', "dna,k=31,scaled=100")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    # check fasta files are present
    fa_files = os.listdir(out_dir)
    print(fa_files)
    assert set(fa_files) == set(['GCA_000175535.1_genomic_first50kb.urlsketch.fna.gz', 'GCA_000175535.1_genomic_second50kb.urlsketch.fna.gz'])

     # Compare the contents of the generated FASTA files to the expected ones
    for generated_file, expected_file in [
        ('GCA_000175535.1_genomic_first50kb.urlsketch.fna.gz', first50kb),
        ('GCA_000175535.1_genomic_second50kb.urlsketch.fna.gz', second50kb)
    ]:
        generated_path = os.path.join(out_dir, generated_file)

        # Read the records from both files using screed
        gen_records = set((record.name, record.sequence) for record in screed.open(generated_path))
        exp_records = set((record.name, record.sequence) for record in screed.open(expected_file))

        # Assert that the records are identical
        assert gen_records == exp_records

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 2
    for sig in sigs:
        ident = sig.name.split(' ')[0]
        assert ident in ["GCA_000175535.1_first50kb", "GCA_000175535.1_second50kb"]
        print(ident)
        if ident == "GCA_000175535.1_first50kb":
            assert sig.md5sum() == ss1.md5sum()
        if ident == "GCA_000175535.1_second50kb":
            assert sig.md5sum() == ss2.md5sum()
    assert os.path.exists(failed)


def test_urlsketch_with_range_improper_range_1(runtmp, capfd):
    acc_csv = get_test_data('acc-url-range.csv')
    acc_mod = runtmp.output("acc-url-range-mod.csv")
    subseqs = get_test_data('subseqs.zip')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')

    # Modify the range in the acc_csv file
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if row['accession'] == 'GCA_000175535.1_second50kb':
                row['range'] = '100000-10000000'
            writer.writerow(row)

    # open subseq sigs
    idx = sourmash.load_file_as_index(subseqs)
    siglist = list(idx.signatures())
    ss1 = siglist[0]

    runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty
    assert os.path.exists(failed)
    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Invalid range: start=100000, end=10000000, sequence length=1088736" in captured.err

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())

    assert len(sigs) == 1
    for sig in sigs:
        ident = sig.name.split(' ')[0]
        assert ident == "GCA_000175535.1_first50kb"
        assert sig.md5sum() == ss1.md5sum()

    with open(failed, 'r') as failF:
        header = next(failF).strip()
        print(header)
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "GCA_000175535.1_second50kb"
            assert name == "GCA_000175535.1_second50kb"
            assert moltype == "DNA"
            assert md5sum == "b9fb20c51f0552b87db5d44d5d4566;a1a8f1c6dc56999c73fe298871c963d1"
            assert download_filename == "both.urlsketch.fna.gz"
            assert url ==  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/961/135/GCA_000961135.2_ASM96113v2/GCA_000961135.2_ASM96113v2_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == ""


def test_urlsketch_with_range_improper_range_2(runtmp, capfd):
    acc_csv = get_test_data('acc-url-range.csv')
    acc_mod = runtmp.output("acc-url-range-mod.csv")
    subseqs = get_test_data('subseqs.zip')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')

    # Modify the range in the acc_csv file
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            if row['accession'] == 'GCA_000175535.1_second50kb':
                row['range'] = '-1-10000000'
            writer.writerow(row)

    # open subseq sigs
    idx = sourmash.load_file_as_index(subseqs)
    siglist = list(idx.signatures())
    ss1 = siglist[0]

    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    captured = capfd.readouterr()
    print(captured.err)
    assert "Error: Invalid range format: -1-10000000" in captured.err


def test_urlsketch_merged_ranged(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    acc_mod = runtmp.output('acc-merged-md5sums-ranges.csv')
    subseqs = get_test_data('subseqs.zip')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')
    sketch_out = runtmp.output('sketch-subseqs.zip')
    merged_out = runtmp.output('merged-subseqs.zip')
    f1 = get_test_data("GCA_000175535.1_ASM17553v1_genomic.1-50000.fna.gz")
    f2 = get_test_data("GCA_000175535.1_ASM17553v1_genomic.50000-100000.fna.gz")

    # Modify the acc_csv file to add range values
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row['url'] = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            row['md5sum'] = "a1a8f1c6dc56999c73fe298871c963d1;a1a8f1c6dc56999c73fe298871c963d1"
            row['range'] = '1-50000;50000-100000'
            writer.writerow(row)
            print(row)

    # sketch subseq files
    runtmp.sourmash('sketch', "dna", f1, f2, '--name',
                    'both name', '-o', sketch_out,
                    '-p', "dna,k=31,scaled=100")

    idx = sourmash.load_file_as_index(sketch_out)
    sigs1 = list(idx.signatures())
    assert len(sigs1) == 1
    sketchsig = sigs1[0]

    # merge subset sketches
    runtmp.sourmash('sig', "merge", subseqs,'--set-name',
                    'both name', '-o', merged_out)
    idx = sourmash.load_file_as_index(merged_out)
    sigs = list(idx.signatures())
    assert len(sigs) == 1
    mergesig = sigs[0]

    # # run urlsketch
    runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    assert os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    idx = sourmash.load_file_as_index(output)
    sigs = list(idx.signatures())
    assert len(sigs) == 1
    sig = sigs[0]
    assert sig.name == "both name"
    print(sig.md5sum())
    assert sig.md5sum() == sketchsig.md5sum() == mergesig.md5sum() == "5feeed4c8a75c8b3fe67af1270fa92c4"


def test_urlsketch_merged_ranged_md5sum_fail_no_checksum_file(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    acc_mod = runtmp.output('acc-merged-md5sums-ranges.csv')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')

    # Modify the acc_csv file to add range values
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row['url'] = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            row['md5sum'] = "a1a8f1c6dc56999c73fe298871c963d1;b2" # second md5sum is incorrect
            row['range'] = '1-50000;50000-100000'
            writer.writerow(row)
            print(row)

    # # run urlsketch
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    assert not os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    with open(failed, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert md5sum == "a1a8f1c6dc56999c73fe298871c963d1;b2"
            assert download_filename == "both.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == "1-50000;50000-100000"


def test_urlsketch_merged_ranged_md5sum_fail_with_checksum_file(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    acc_mod = runtmp.output('acc-merged-md5sums-ranges.csv')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')
    ch_fail = runtmp.output('ch_failed.csv')

    # Modify the acc_csv file to add range values
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            row['url'] = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            row['md5sum'] = "a1a8f1c6dc56999c73fe298871c963d1;b2" # second md5sum is incorrect
            row['range'] = '1-50000;50000-100000'
            writer.writerow(row)
            print(row)

    # # run urlsketch
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1', '--checksum-fail', ch_fail,
                    '--param-str', "dna,k=31,scaled=100")

    assert not os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    # since this is a merged dataset, we write both checksum fail and regular fail.
    with open(failed, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert md5sum == "a1a8f1c6dc56999c73fe298871c963d1;b2"
            assert download_filename == "both.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == "1-50000;50000-100000"

    assert os.path.exists(ch_fail)
    with open(ch_fail, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum_url,download_filename,url,expected_md5sum,reason"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum_url, download_filename, url, expected_md5, reason= line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert md5sum_url == ""
            assert download_filename == "both.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert expected_md5 == "b2"
            assert reason == "MD5 hash does not match. Expected: 'b2'; Found: 'a1a8f1c6dc56999c73fe298871c963d1'"


def test_urlsketch_merged_ranged_fail(runtmp):
    acc_csv = get_test_data('acc-merged-md5sums.csv')
    acc_mod = runtmp.output('acc-merged-md5sums-ranges.csv')
    output = runtmp.output('range.zip')
    failed = runtmp.output('failed.csv')

    # Modify the acc_csv file to add range values
    with open(acc_csv, 'r') as infile, open(acc_mod, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader:
            # first url is incorrect
            row['url'] = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            row['md5sum'] = "a1a8f1c6dc56999c73fe298871c963d1;"
            row['range'] = '1-50000;50000-100000'
            writer.writerow(row)
            print(row)

    # # run urlsketch
    with pytest.raises(utils.SourmashCommandFailed):
        runtmp.sourmash('scripts', 'urlsketch', acc_mod, '-o', output,
                    '--failed', failed, '-r', '1',
                    '--param-str', "dna,k=31,scaled=100")

    assert not os.path.exists(output)
    assert not runtmp.last_result.out # stdout should be empty

    with open(failed, 'r') as failF:
        header = next(failF).strip()
        assert header == "accession,name,moltype,md5sum,download_filename,url,range"
        for line in failF:
            print(line)
            acc, name, moltype, md5sum, download_filename, url, range = line.strip().split(',')
            assert acc == "both"
            assert name == "both name"
            assert moltype == "DNA"
            assert md5sum == "a1a8f1c6dc56999c73fe298871c963d1;"
            assert download_filename == "both.urlsketch.fna.gz"
            assert url == "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1;https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/175/535/GCA_000175535.1_ASM17553v1/GCA_000175535.1_ASM17553v1_genomic.fna.gz"
            assert range == "1-50000;50000-100000"
