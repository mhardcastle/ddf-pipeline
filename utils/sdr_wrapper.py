"""LoTSS-DR3 wrapper with the historical SDR calling interface.

get_status returns DUL (online), OFL (offline), or STG (restoring).
stage returns 'Online' or 'Staged'; download_and_stage requests restoration,
polls until online, then downloads. Requires requests and ddf-pipeline's
download_file.py helper.

SURF's Request button GETs its web download URL to initiate restoration.
Archived readiness is read from the same server-rendered file table; the
upload status 'completed' is NOT evidence that a file is online.

Environment: SDR_URL (records API root), SDR_TOKEN (new SURF bearer token),
SDR_POLL_INTERVAL (seconds, default 30), SDR_STAGE_TIMEOUT (seconds,
default 0 = wait indefinitely, matching the old wrapper).
"""

from __future__ import print_function

import csv
import io
import os
import random
import re
import uuid
from html.parser import HTMLParser
from time import monotonic, sleep
from urllib.parse import quote, unquote, urljoin, urlsplit

import requests
from download_file import download_file


def request_with_retry(url, rfunction=requests.get, **kwargs):
    """Retry transient failures, not authentication or certificate failures."""
    kwargs.setdefault('timeout', 60)
    kwargs.setdefault('verify', True)
    for attempt in range(6):
        try:
            response = rfunction(url, **kwargs)
        except requests.exceptions.SSLError:
            raise
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            if attempt == 5:
                raise
        else:
            if response.status_code not in (429, 500, 502, 503, 504) or attempt == 5:
                return response
            response.close()
        sleep(min(30, 2 ** attempt) + random.uniform(0, 2))


class _BearerAuth(requests.auth.AuthBase):
    """Never send a repository token to ASTRON or another origin."""
    def __init__(self, token, origin):
        self.token, self.origin = token, origin

    def __call__(self, request):
        parsed = urlsplit(request.url)
        if self.token and (parsed.scheme, parsed.netloc) == self.origin:
            request.headers['Authorization'] = 'Bearer ' + self.token
        return request


class _FileTableParser(HTMLParser):
    """Read current labels and restore links from SURF's own file table."""
    def __init__(self):
        HTMLParser.__init__(self, convert_charrefs=True)
        self.files = {}
        self._section = 0
        self._row = None
        self._name_cell = self._name_p = self._label = False

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        classes = (attrs.get('class') or '').split()
        if tag == 'section':
            if self._section:
                self._section += 1
            elif attrs.get('id') == 'record-files':
                self._section = 1
        if not self._section:
            return
        if tag == 'tr':
            self._row = None if 'accordion-row' in classes else {
                'name': [], 'label': [], 'main': False, 'key': None,
                'request_url': None, 'download_url': None}
        if self._row is None:
            return
        if tag == 'td':
            self._name_cell = 'collapsible-name-cell' in classes
            self._row['main'] |= self._name_cell
        if tag == 'p' and self._name_cell:
            self._name_p = True
        if tag == 'div' and 'label' in classes:
            self._label = True
        if tag == 'a':
            disabled = 'disabled' in classes or attrs.get('aria-disabled') == 'true'
            if 'request-file-button' in classes:
                self._row['key'] = attrs.get('data-file-key')
                if not disabled:
                    self._row['request_url'] = attrs.get('data-request-url')
            elif not disabled and '/files/' in (attrs.get('href') or ''):
                self._row['download_url'] = attrs['href']

    def handle_data(self, text):
        if self._row is not None:
            if self._name_p:
                self._row['name'].append(text)
            if self._label:
                self._row['label'].append(text)

    def handle_endtag(self, tag):
        if tag == 'p':
            self._name_p = False
        if tag == 'td':
            self._name_cell = False
        if tag == 'div':
            self._label = False
        if tag == 'tr' and self._row is not None:
            row = self._row
            name = ''.join(row['name']).strip()
            label = ''.join(row['label']).strip().lower()
            status = None
            if label == 'online' and row['download_url']:
                status = 'DUL'
            elif label == 'external' and row['download_url']:
                status = 'EXT'  # Resolve the linked file; this is not an online label.
            elif label == 'offline' and row['request_url']:
                status = 'OFL'
            elif any(w in label for w in ('restor', 'stag', 'request', 'pending', 'retriev')):
                status = 'STG'
            if row['main'] and name and status:
                if row['key'] and row['key'] != name:
                    raise RuntimeError('Inconsistent filename in SURF file table')
                self.files[name] = {
                    'status': status, 'request_url': row['request_url'],
                    'download_url': row['download_url']}
            self._row = None
        if tag == 'section' and self._section:
            self._section -= 1


class SDR(object):
    TAP_URL = 'https://vo.astron.nl/tap/sync'

    def __init__(self, url='https://repository.surf.nl/api/records/', target='.'):
        self.url = os.environ.get('SDR_URL', url).rstrip('/') + '/'
        parsed = urlsplit(self.url)
        if parsed.scheme != 'https' or not parsed.netloc:
            raise ValueError('SDR_URL must be an HTTPS records API URL')
        if '/api/objects/' in parsed.path:
            raise ValueError('SDR_URL points to the old API. Unset it or use '
                             'https://repository.surf.nl/api/records/')
        self.token = os.environ.get('SDR_TOKEN')
        self.tokenstr = ''  # Compatibility attribute; no DR2 share-token in DR3.
        self.target = target
        self.session = requests.Session()
        self._auth = _BearerAuth(self.token, (parsed.scheme, parsed.netloc))
        self.session.auth = self._auth
        self._records = {}
        self.poll_interval = float(os.environ.get('SDR_POLL_INTERVAL', '30'))
        self.stage_timeout = float(os.environ.get('SDR_STAGE_TIMEOUT', '0'))
        if self.poll_interval < 0 or self.stage_timeout < 0:
            raise ValueError('Polling interval and timeout must not be negative')

    def fc(self, field):
        """Retain old field conversion. Both P213+73 and P213-73 are accepted."""
        return field.replace('+', '-')

    @staticmethod
    def _check_response(response):
        if response.status_code in (401, 403):
            raise RuntimeError('SURF denied access (HTTP %s). Public DR3 data '
                               'normally need no token; otherwise use a new '
                               'SURF bearer token, not a DR2 share token.' %
                               response.status_code)
        response.raise_for_status()
        if re.search(r'/(login|signin|oauth|saml)(/|$)', urlsplit(response.url).path):
            raise RuntimeError('The repository redirected to a sign-in page')

    def _json(self, url):
        response = request_with_retry(url, rfunction=self.session.get,
                                      headers={'Cache-Control': 'no-cache'})
        try:
            self._check_response(response)
            return response.json()
        finally:
            response.close()

    def _get_record(self, field):
        if field in self._records:
            return self._records[field]
        names = [field]
        if re.fullmatch(r'P\d+-\d+', field):
            names.append(field.replace('-', '+', 1))
        values = ','.join("'%s'" % name.replace("'", "''") for name in names)
        query = ('SELECT pointing_id, record_id FROM lotss_dr3.pointings '
                 'WHERE pointing_id IN (%s)' % values)
        # Do not send the SURF token with the ASTRON lookup.
        response = request_with_retry(self.TAP_URL, params={
            'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': query})
        try:
            response.raise_for_status()
            reader = csv.DictReader(io.StringIO(response.text.lstrip('\ufeff')))
            if not reader.fieldnames or not {'pointing_id', 'record_id'}.issubset(reader.fieldnames):
                raise RuntimeError('Unexpected ASTRON TAP response; expected CSV')
            rows = list(reader)
        finally:
            response.close()
        exact = [row for row in rows if row['pointing_id'] == field]
        if exact:
            rows = exact
        if len(rows) != 1 or not rows[0]['record_id'].strip():
            raise RuntimeError('DR3 pointing not found or ambiguous: %s' % field)
        record = self._json(self.url + quote(rows[0]['record_id'].strip(), safe=''))
        self._records[field] = record  # Identity only, not availability.
        return record

    def _page_files(self, record):
        page_url = record['links']['self_html']
        response = request_with_retry(page_url, rfunction=self.session.get,
                                      headers={'Cache-Control': 'no-cache'})
        try:
            self._check_response(response)
            parser = _FileTableParser()
            parser.feed(response.text)
        finally:
            response.close()
        for action in parser.files.values():
            for key in ('request_url', 'download_url'):
                if action[key]:
                    action[key] = urljoin(page_url, action[key])
        return parser.files

    def _entries(self, field):
        return self._record_entries(self._get_record(field))

    def _record_entries(self, record):
        listing = self._json(record['links']['files'])
        raw = listing.get('entries')
        if isinstance(raw, dict):
            return {key: dict(value, key=value.get('key', key))
                    for key, value in raw.items()}
        if isinstance(raw, list):
            return {entry['key']: entry for entry in raw}
        raise RuntimeError('Unexpected SURF file-list response')

    def _external_file(self, entry, records):
        """Resolve published external links without GETting/staging the linked data.

        Some DR3 files point via EPIC handles to archived files in DR2 records.
        Stop at the target SURF file URL and read its live metadata/file table.
        The per-snapshot cache must not retain availability across polling calls.
        """
        url = urljoin(self.url, entry['links']['content'])
        visited = set()
        origin = urlsplit(self.url)
        for _ in range(10):
            if url in visited:
                break
            visited.add(url)
            response = request_with_retry(url, rfunction=self.session.head,
                                          allow_redirects=False)
            try:
                self._check_response(response)
                if response.status_code in (301, 302, 303, 307, 308):
                    location = response.headers.get('Location')
                    if not location:
                        raise RuntimeError('External file redirect has no Location')
                    url = urljoin(response.url, location)
                elif (response.status_code == 200 and
                      response.headers.get('Content-Length') == str(entry['size']) and
                      response.headers.get('Content-Type', '').split(';')[0].strip().lower()
                      not in ('', 'text/html', 'application/json')):
                    return entry, {'status': 'DUL', 'request_url': None,
                                   'download_url': url}
                else:
                    raise RuntimeError('Cannot determine external availability for %s '
                                       '(HTTP %s)' % (entry['key'], response.status_code))
            finally:
                response.close()
            target = urlsplit(url)
            match = re.search(r'/(?:api/)?records/([^/]+)/files/([^/]+)(?:/content)?/?$', target.path)
            if (target.scheme, target.netloc) == (origin.scheme, origin.netloc) and match:
                record_id, filename = (unquote(value) for value in match.groups())
                if record_id not in records:
                    record = self._json(self.url + quote(record_id, safe=''))
                    records[record_id] = (self._record_entries(record), None, record)
                entries, actions, record = records[record_id]
                if filename not in entries:
                    raise RuntimeError('External file missing from linked SURF record: %s' % filename)
                entry = entries[filename]
                if str(entry.get('storage_status', '')).upper() == 'STANDARD':
                    return entry, {'status': 'DUL', 'request_url': None,
                                   'download_url': entry['links']['content']}
                if actions is None:
                    actions = self._page_files(record)
                    records[record_id] = (entries, actions, record)
                action = actions.get(filename)
                if action and action['status'] != 'EXT':
                    return entry, action
                if not action:
                    raise RuntimeError('Cannot determine linked SURF availability for %s' % filename)
                url = urljoin(self.url, entry['links']['content'])
        raise RuntimeError('External SURF file link contains a redirect loop: %s' % entry['key'])

    def _snapshot(self, field):
        entries = self._entries(field)
        statuses, actions = {}, {}
        for name, entry in entries.items():
            if str(entry.get('storage_status', '')).upper() == 'STANDARD':
                statuses[name] = 'DUL'
        if len(statuses) != len(entries):
            # The restore boolean alone does not establish readiness. Use the
            # live UI's Online/Offline label, including restored GLACIER files.
            actions = self._page_files(self._get_record(field))
            external_records = {}
            for name in entries:
                if name in actions:
                    if actions[name]['status'] == 'EXT':
                        entry, action = self._external_file(entries[name], external_records)
                        entries[name] = dict(entry, key=name)
                        actions[name] = action
                    statuses[name] = actions[name]['status']
                elif name not in statuses:
                    raise RuntimeError('Cannot determine storage availability for %s; '
                                       'the SURF file table may have changed.' % name)
        return entries, statuses, actions

    def get_files(self, field):
        return {name: i for i, name in enumerate(self._entries(field))}

    def get_status(self, field):
        return self._snapshot(field)[1]

    def stage(self, field, filename):
        entries, statuses, actions = self._snapshot(field)
        if filename not in entries:
            raise RuntimeError('Staging nonexistent file: %s' % filename)
        if statuses[filename] == 'DUL':
            return 'Online'
        if statuses[filename] == 'STG':
            return 'Staged'
        request_url = actions.get(filename, {}).get('request_url')
        if not request_url:
            raise RuntimeError('No restore-request link for %s' % filename)
        # Match SURF's Request button. Do not consume a huge file body if the
        # restore completes between the status check and this GET.
        response = request_with_retry(request_url, rfunction=self.session.get, stream=True)
        try:
            self._check_response(response)
        finally:
            response.close()
        return 'Staged'

    def _wait_online(self, field, filenames):
        started = monotonic()
        print('Waiting for files to be online:')
        while True:
            statuses = self.get_status(field)
            missing = [name for name in filenames if name not in statuses]
            if missing:
                raise RuntimeError('Files disappeared from repository: ' + ', '.join(missing))
            ready = sum(statuses[name] == 'DUL' for name in filenames)
            print('%i/%i... ' % (ready, len(filenames)), end='', flush=True)
            if ready == len(filenames):
                print()
                return
            elapsed = monotonic() - started
            if self.stage_timeout and elapsed >= self.stage_timeout:
                raise RuntimeError('Timed out waiting for SURF to restore files')
            delay = self.poll_interval
            if self.stage_timeout:
                delay = min(delay, max(0, self.stage_timeout - elapsed))
            sleep(delay)

    def download(self, field, filename, progress_bar=False):
        entries, statuses, _ = self._snapshot(field)
        if filename not in entries:
            raise RuntimeError('File not found: %s' % filename)
        if statuses[filename] != 'DUL':
            raise RuntimeError('File not online: %s' % filename)
        entry = entries[filename]
        url = entry.get('links', {}).get('content')
        if not url:
            raise RuntimeError('No content link for %s' % filename)
        url = urljoin(self.url, url)
        root = os.path.realpath(self.target)
        destination = os.path.realpath(os.path.join(root, filename))
        if os.path.commonpath([root, destination]) != root or destination == root:
            raise RuntimeError('Unsafe filename: %s' % filename)
        os.makedirs(os.path.dirname(destination), exist_ok=True)
        partial = destination + '.part'
        if os.path.islink(partial):
            raise RuntimeError('Partial-download path must not be a symlink')
        expected_size = int(entry['size'])
        if os.path.isfile(destination) and os.path.getsize(destination) == expected_size:
            print('Already downloaded:', filename)
            return
        partial_size = os.path.getsize(partial) if os.path.exists(partial) else 0
        if partial_size > expected_size:
            raise RuntimeError('Existing partial file is too large: %s' % partial)
        if partial_size != expected_size or not os.path.exists(partial):
            download_file(url, partial, catch_codes=(429, 500, 502, 503, 504),
                          retry_partial=True, force_retry_partial=partial_size > 0,
                          progress_bar=progress_bar, verify=True,
                          auth=self._auth, maxsleep=10)
        if not os.path.isfile(partial):
            raise RuntimeError('Downloader did not create the expected file: %s' % partial)
        if os.path.getsize(partial) != expected_size:
            # A successful HTTP 200 can still be a small restoration/error
            # response after a storage race. Do not resume that as a data prefix.
            rejected = partial + '.invalid-' + uuid.uuid4().hex
            os.replace(partial, rejected)
            raise RuntimeError('Downloaded size differs from SURF metadata for %s; '
                               'unexpected response preserved at %s' % (filename, rejected))
        os.replace(partial, destination)

    def download_and_stage(self, field, filenames, progress_bar=False):
        filenames = list(dict.fromkeys(filenames))
        if not filenames:
            return
        statuses = self.get_status(field)
        missing = [name for name in filenames if name not in statuses]
        if missing:
            raise RuntimeError('File not found: ' + ', '.join(missing))
        for filename in filenames:
            if statuses[filename] == 'OFL':
                print('Requesting restore:', filename)
                self.stage(field, filename)
        self._wait_online(field, filenames)
        for filename in filenames:
            # Large downloads may outlast the restoration window of later files.
            current = self.get_status(field)
            if filename not in current:
                raise RuntimeError('File disappeared from repository: %s' % filename)
            if current[filename] != 'DUL':
                if current[filename] == 'OFL':
                    self.stage(field, filename)
                self._wait_online(field, [filename])
            print('Downloading', filename)
            self.download(field, filename, progress_bar=progress_bar)

