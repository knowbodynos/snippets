import boto3
import botocore

def pull(bucket_name, s3_keys):
    """
    Download files from s3 bucket

    :param str bucket_name: name of s3 bucket
    :param dict s3_keys: dictionary with s3 keys and the corresponding local paths to download them to
    """
    s3 = boto3.resource('s3')
    try:
        for s3_key, local_path in s3_keys.items():
            s3.Bucket(bucket_name).download_file(s3_key, local_path)
    except botocore.exceptions.ClientError as e:
        if e.response['Error']['Code'] == "404":
            print("The object does not exist.")
        else:
            raise