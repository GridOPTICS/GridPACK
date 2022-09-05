#!/usr/bin/env python

import sys

def github_api_json(address):
    import json
    if (sys.version_info > (3, 0)):
        import urllib.request
        with urllib.request.urlopen(address) as url:
            data = json.loads(url.read().decode())
    else:
        import urllib
        response = urllib.urlopen(address)
        data = json.loads(response.read())
    return data

if len(sys.argv) != 3:
    print("usage: github_download_counter.py <organization> <repo>")
    sys.exit(1)

org = sys.argv[1]
repo = sys.argv[2]

address = 'https://api.github.com/repos/{}/{}/releases'.format(org,repo)

print('Download counts for {}/{}'.format(org,repo))

releases = github_api_json(address)
if not releases or 'message' in releases:
    print("failed to download github api release JSON")
    sys.exit(1)

for release in releases:
    print(release['tag_name'])
    for asset in release['assets']:
        print('    {}: {}'.format(asset['name'], asset['download_count']))
