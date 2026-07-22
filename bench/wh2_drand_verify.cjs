#!/usr/bin/env node
"use strict";

// Entropy/root acquisition always asks every pinned origin for one
// caller-selected round.  The separate latest-bound mode only establishes a
// verified pre-seal clock/liveness lower bound; its output is never entropy.

const crypto = require("crypto");
const fs = require("fs");
const path = require("path");

const CHAIN_HASH =
  "52db9ba70e0cc0f6eaf7803dd07447a1f5477735fd3f661792ba94600c84e971";
const PUBLIC_KEY =
  "83cf0f2896adee7eb8b5f01fcad3912212c437e0073e911fb90022d3e760183c8" +
  "c4b450b6a0a6c3ac6a5776a2d1064510d1fec758c921cc22b0e17e63aaf4bcb5" +
  "ed66304de9cf809bd274ca73bab4af5a6e9c76a4bc09e76eae8991ef5ece45a";
const PINNED_INFO = Object.freeze({
  public_key: PUBLIC_KEY,
  period: 3,
  genesis_time: 1692803367,
  hash: CHAIN_HASH,
  groupHash:
    "f477d5c89f21a17c863a7f937c6a6d15859414d2be09cd448d4279af331c5d3e",
  schemeID: "bls-unchained-g1-rfc9380",
  metadata: Object.freeze({beaconID: "quicknet"}),
});
const ORIGINS = Object.freeze([
  "https://api.drand.sh",
  "https://api2.drand.sh",
  "https://api3.drand.sh",
  "https://drand.cloudflare.com",
]);
const TIMEOUT_MS = 8000;
const RESPONSE_LIMIT_BYTES = 4096;
const SELFTEST = Object.freeze({
  round: 1,
  randomness:
    "1466a6cd24e327188770752f6134001c64d6efcc590ccc26b721611ad96f165a",
  signature:
    "b55e7cb2d5c613ee0b2e28d6750aabbb78c39dcc96bd9d38c2c2e12198df955" +
    "71de8e8e402a0cc48871c7089a2b3af4b",
});
const SELFTEST_CANONICAL_SHA256 =
  "b03d725b13e383cb643e02e448bf663e70e0f1b9dc85c33b0e34cefe151efeff";

function fail(message) {
  throw new Error(message);
}

function parseRound(text) {
  if (!/^[1-9][0-9]*$/.test(text)) {
    fail("round must be a canonical positive decimal integer");
  }
  const value = Number(text);
  if (!Number.isSafeInteger(value)) {
    fail("round does not fit a JavaScript safe integer");
  }
  return value;
}

function canonicalBeacon(beacon, expectedRound, context) {
  if (beacon === null || typeof beacon !== "object" ||
      Array.isArray(beacon) || Object.getPrototypeOf(beacon) !== Object.prototype ||
      !Number.isSafeInteger(beacon.round) ||
      beacon.round !== expectedRound ||
      typeof beacon.randomness !== "string" ||
      !/^[0-9a-fA-F]{64}$/.test(beacon.randomness) ||
      typeof beacon.signature !== "string" ||
      !/^[0-9a-fA-F]{96}$/.test(beacon.signature) ||
      Object.prototype.hasOwnProperty.call(beacon, "previous_signature")) {
    fail(`${context}: invalid Quicknet beacon semantics`);
  }
  const normalized = {
    round: beacon.round,
    randomness: beacon.randomness.toLowerCase(),
    signature: beacon.signature.toLowerCase(),
  };
  const bytes = Buffer.from(JSON.stringify(normalized) + "\n", "ascii");
  return {normalized, bytes};
}

function pinnedChain(baseUrl, info = PINNED_INFO) {
  return {baseUrl, info: async () => info};
}

function clientOptions(drand) {
  return {
    ...drand.defaultChainOptions,
    disableBeaconVerification: false,
    noCache: true,
    chainVerificationParams: {chainHash: CHAIN_HASH, publicKey: PUBLIC_KEY},
  };
}

function classifyLatestObservations(observations) {
  const groups = new Map();
  const canonicalKeysByRound = new Map();
  for (const item of observations.filter((item) => item.ok)) {
    const key = JSON.stringify(item.beacon);
    const members = groups.get(key) || [];
    members.push(item.origin);
    groups.set(key, members);
    const roundKeys = canonicalKeysByRound.get(item.beacon.round) || new Set();
    roundKeys.add(key);
    canonicalKeysByRound.set(item.beacon.round, roundKeys);
  }
  // Adjacent latest rounds are normal at pulse boundaries.  Two different
  // BLS-verified canonical values for the same exact round are not.
  if ([...canonicalKeysByRound.values()].some((keys) => keys.size > 1)) {
    return {status: "PERMANENT_VERIFIED_DISAGREEMENT"};
  }
  const consensus = [...groups.entries()].filter((entry) => entry[1].length >= 3);
  if (consensus.length === 1) {
    return {
      status: "QUORUM",
      canonicalKey: consensus[0][0],
      consensusOrigins: consensus[0][1],
    };
  }
  if (consensus.length > 1) {
    return {status: "PERMANENT_VERIFIED_DISAGREEMENT"};
  }
  return {status: "TEMPORARY_NO_QUORUM"};
}

function latestGroupingKats() {
  const newer = {...SELFTEST, round: SELFTEST.round + 1};
  const sameRoundConflict = {
    ...SELFTEST,
    randomness: `${SELFTEST.randomness[0] === "0" ? "1" : "0"}` +
      SELFTEST.randomness.slice(1),
  };
  const observation = (origin, beacon) => ({origin, ok: true, beacon});
  const threePlusOne = classifyLatestObservations([
    observation(ORIGINS[0], SELFTEST),
    ...ORIGINS.slice(1).map((origin) => observation(origin, newer)),
  ]);
  const twoPlusTwo = classifyLatestObservations([
    ...ORIGINS.slice(0, 2).map((origin) => observation(origin, SELFTEST)),
    ...ORIGINS.slice(2).map((origin) => observation(origin, newer)),
  ]);
  const impossibleSameRound = classifyLatestObservations([
    ...ORIGINS.slice(0, 2).map((origin) => observation(origin, SELFTEST)),
    ...ORIGINS.slice(2).map(
      (origin) => observation(origin, sameRoundConflict)),
  ]);
  if (threePlusOne.status !== "QUORUM" ||
      JSON.parse(threePlusOne.canonicalKey).round !== newer.round ||
      threePlusOne.consensusOrigins.length !== 3 ||
      twoPlusTwo.status !== "TEMPORARY_NO_QUORUM" ||
      impossibleSameRound.status !== "PERMANENT_VERIFIED_DISAGREEMENT") {
    fail("latest-bound grouping KAT failed");
  }
  return [
    "adjacent-round-3-plus-1-quorum",
    "adjacent-round-2-plus-2-temporary",
    "same-round-verified-split-permanent",
  ];
}

async function expectReject(operation, context) {
  const savedConsoleError = console.error;
  console.error = () => {};
  try {
    await operation();
  } catch (_) {
    return;
  } finally {
    console.error = savedConsoleError;
  }
  fail(`offline mutation KAT unexpectedly accepted ${context}`);
}

async function offlineSelftest(drand) {
  if (typeof drand.fetchBeacon !== "function") {
    fail("pinned drand-client lacks required public exports");
  }
  const fake = (beacon, info = PINNED_INFO) => ({
    options: clientOptions(drand),
    get: async () => beacon,
    latest: async () => fail("latest is forbidden"),
    chain: () => pinnedChain("offline://quicknet", info),
  });
  const accepted = await drand.fetchBeacon(fake(SELFTEST), SELFTEST.round);
  const canonical = canonicalBeacon(accepted, SELFTEST.round, "offline KAT");
  const canonicalSha = crypto.createHash("sha256").update(canonical.bytes).digest("hex");
  if (canonicalSha !== SELFTEST_CANONICAL_SHA256) {
    fail("offline canonical-beacon SHA256 golden mismatch");
  }
  const mutate = (text) => `${text[0] === "0" ? "1" : "0"}${text.slice(1)}`;
  await expectReject(
    () => drand.fetchBeacon(fake({...SELFTEST, randomness: mutate(SELFTEST.randomness)}), 1),
    "randomness",
  );
  await expectReject(
    () => drand.fetchBeacon(fake({...SELFTEST, signature: mutate(SELFTEST.signature)}), 1),
    "signature",
  );
  await expectReject(() => drand.fetchBeacon(fake(SELFTEST), 2), "expected round");
  await expectReject(
    () => drand.fetchBeacon(fake(SELFTEST, {...PINNED_INFO, public_key: mutate(PUBLIC_KEY)}), 1),
    "public key",
  );
  await expectReject(
    () => drand.fetchBeacon(fake(SELFTEST, {...PINNED_INFO, schemeID: "invalid"}), 1),
    "scheme",
  );
  const groupingKats = latestGroupingKats();
  return {
    schema: "wirehair.wh2.drand_quicknet_offline_kat.v1",
    beacon: canonical.normalized,
    canonical_sha256: canonicalSha,
    mutation_rejections: ["randomness", "signature", "round", "key", "scheme"],
    latest_bound_grouping_kats: groupingKats,
  };
}

async function offlineVerify(drand) {
  const raw = fs.readFileSync(0);
  if (raw.length > RESPONSE_LIMIT_BYTES) {
    fail("offline beacon input exceeds 4096 bytes");
  }
  let beacon;
  try {
    beacon = JSON.parse(raw.toString("utf8"));
  } catch (_) {
    fail("offline beacon input is malformed JSON");
  }
  const canonicalInput = canonicalBeacon(beacon, beacon.round, "offline verify");
  const client = {
    options: clientOptions(drand),
    get: async (round) => {
      if (round !== beacon.round) fail("offline verifier requested another round");
      return beacon;
    },
    latest: async () => fail("latest is forbidden"),
    chain: () => pinnedChain("offline://quicknet"),
  };
  const verified = await drand.fetchBeacon(client, beacon.round);
  const canonical = canonicalBeacon(verified, beacon.round, "offline verify");
  if (!canonical.bytes.equals(canonicalInput.bytes)) {
    fail("offline verifier changed canonical beacon bytes");
  }
  return {
    schema: "wirehair.wh2.drand_quicknet_offline_verified.v1",
    beacon: canonical.normalized,
    canonical_sha256:
      crypto.createHash("sha256").update(canonical.bytes).digest("hex"),
  };
}

async function fetchBoundedResponse(nativeFetch, origin, url) {
  const startedMs = Date.now();
  const response = await nativeFetch(url, {
    method: "GET",
    redirect: "error",
    signal: AbortSignal.timeout(TIMEOUT_MS),
    headers: {Accept: "application/json"},
  });
  if (!response.ok || response.status !== 200 || response.redirected ||
      response.url !== url) {
    fail(`${origin}: noncanonical HTTP response`);
  }
  if (!response.body || typeof response.body.getReader !== "function") {
    fail(`${origin}: response has no bounded web stream`);
  }
  const reader = response.body.getReader();
  const chunks = [];
  let size = 0;
  for (;;) {
    const {done, value} = await reader.read();
    if (done) break;
    size += value.byteLength;
    if (size > RESPONSE_LIMIT_BYTES) {
      await reader.cancel();
      fail(`${origin}: response exceeds ${RESPONSE_LIMIT_BYTES} bytes`);
    }
    chunks.push(Buffer.from(value));
  }
  const raw = Buffer.concat(chunks, size);
  return {
    raw,
    transport: {
      fetch_started_ms: startedMs,
      fetch_completed_ms: Date.now(),
      raw_response_sha256:
        crypto.createHash("sha256").update(raw).digest("hex"),
    },
  };
}

async function fetchRound(drand, round) {
  const nativeFetch = globalThis.fetch;
  if (typeof nativeFetch !== "function") {
    fail("pinned Node runtime does not provide native fetch");
  }
  const observations = await Promise.all(ORIGINS.map(async (origin) => {
    let transport = {};
    try {
      const chain = pinnedChain(`${origin}/${CHAIN_HASH}`);
      const client = {
        options: clientOptions(drand),
        chain: () => chain,
        latest: async () => fail("latest beacon requests are forbidden"),
        get: async (requestedRound) => {
          if (requestedRound !== round) {
            fail(`${origin}: verifier requested an uncommitted round`);
          }
          const url = `${origin}/${CHAIN_HASH}/public/${round}`;
          const bounded = await fetchBoundedResponse(nativeFetch, origin, url);
          const raw = bounded.raw;
          transport = bounded.transport;
          let parsed;
          try {
            parsed = JSON.parse(raw.toString("utf8"));
          } catch (_) {
            fail(`${origin}: malformed beacon JSON`);
          }
          return parsed;
        },
      };
      const verified = await drand.fetchBeacon(client, round);
      const canonical = canonicalBeacon(verified, round, origin);
      return {
        origin,
        ok: true,
        canonical_sha256:
          crypto.createHash("sha256").update(canonical.bytes).digest("hex"),
        beacon: canonical.normalized,
        ...transport,
      };
    } catch (error) {
      return {origin, ok: false, ...transport, error: String(error.message || error)};
    }
  }));
  const classification = classifyLatestObservations(observations);
  if (classification.status === "PERMANENT_VERIFIED_DISAGREEMENT") {
    return {
      schema: "wirehair.wh2.drand_quicknet_wave.v1",
      status: "PERMANENT_VERIFIED_DISAGREEMENT",
      chain_hash: CHAIN_HASH,
      observations,
    };
  }
  if (classification.status === "TEMPORARY_NO_QUORUM") {
    return {
      schema: "wirehair.wh2.drand_quicknet_wave.v1",
      status: "TEMPORARY_NO_QUORUM",
      chain_hash: CHAIN_HASH,
      observations,
    };
  }
  return {
    schema: "wirehair.wh2.drand_quicknet_wave.v1",
    status: "QUORUM",
    chain_hash: CHAIN_HASH,
    consensus_origins: classification.consensusOrigins,
    observations,
    beacon: JSON.parse(classification.canonicalKey),
  };
}

async function fetchLatestBound(drand) {
  const nativeFetch = globalThis.fetch;
  if (typeof nativeFetch !== "function") {
    fail("pinned Node runtime does not provide native fetch");
  }
  // drand-client reports rejected signatures to console.error before it
  // rejects.  Keep all per-origin failures in the structured observation
  // ledger so the wrapper's stderr remains reserved for process-level faults.
  const savedConsoleError = console.error;
  let observations;
  console.error = () => {};
  try {
    observations = await Promise.all(ORIGINS.map(async (origin) => {
      let transport = {};
      try {
        const url = `${origin}/${CHAIN_HASH}/public/latest`;
        const bounded = await fetchBoundedResponse(nativeFetch, origin, url);
        const raw = bounded.raw;
        transport = bounded.transport;
        let parsed;
        try {
          parsed = JSON.parse(raw.toString("utf8"));
        } catch (_) {
          fail(`${origin}: malformed latest-beacon JSON`);
        }
        if (parsed === null || typeof parsed !== "object" ||
            Array.isArray(parsed) ||
            !Number.isSafeInteger(parsed.round) || parsed.round < 1) {
          fail(`${origin}: latest beacon has a noncanonical round`);
        }
        const round = parsed.round;
        const client = {
          options: clientOptions(drand),
          chain: () => pinnedChain(`${origin}/${CHAIN_HASH}`),
          latest: async () => fail("drand-client latest requests are forbidden"),
          get: async (requestedRound) => {
            if (requestedRound !== round) {
              fail(`${origin}: latest-bound verifier requested another round`);
            }
            return parsed;
          },
        };
        const verified = await drand.fetchBeacon(client, round);
        const canonical = canonicalBeacon(verified, round, origin);
        return {
          origin,
          ok: true,
          canonical_sha256:
            crypto.createHash("sha256").update(canonical.bytes).digest("hex"),
          beacon: canonical.normalized,
          ...transport,
        };
      } catch (error) {
        return {
          origin,
          ok: false,
          ...transport,
          error: String(error && error.message ? error.message : error),
        };
      }
    }));
  } finally {
    console.error = savedConsoleError;
  }
  const classification = classifyLatestObservations(observations);
  const base = {
    schema: "wirehair.wh2.drand_quicknet_latest_wave.v1",
    role: "preseal-clock-liveness-lower-bound-only",
    chain_hash: CHAIN_HASH,
  };
  if (classification.status === "PERMANENT_VERIFIED_DISAGREEMENT") {
    return {
      ...base,
      status: "PERMANENT_VERIFIED_DISAGREEMENT",
      observations,
    };
  }
  if (classification.status === "TEMPORARY_NO_QUORUM") {
    return {
      ...base,
      status: "TEMPORARY_NO_QUORUM",
      observations,
    };
  }
  const beacon = JSON.parse(classification.canonicalKey);
  return {
    ...base,
    status: "QUORUM",
    latest_round: beacon.round,
    consensus_origins: classification.consensusOrigins,
    observations,
    beacon,
  };
}

async function main(argv) {
  if (argv.length !== 2) {
    fail("usage: wh2_drand_verify.cjs DRAND_CLIENT_CJS offline-selftest|offline-verify|latest-bound|ROUND");
  }
  const bundle = path.resolve(argv[0]);
  const stat = fs.lstatSync(bundle);
  if (!stat.isFile() || stat.isSymbolicLink()) {
    fail("drand-client bundle is not a regular nonsymlink file");
  }
  const drand = require(bundle);
  let record;
  if (argv[1] === "offline-selftest") {
    record = await offlineSelftest(drand);
  } else if (argv[1] === "offline-verify") {
    record = await offlineVerify(drand);
  } else if (argv[1] === "latest-bound") {
    record = await fetchLatestBound(drand);
  } else {
    record = await fetchRound(drand, parseRound(argv[1]));
  }
  process.stdout.write(JSON.stringify(record) + "\n");
  if (record.status === "TEMPORARY_NO_QUORUM") process.exitCode = 2;
  if (record.status === "PERMANENT_VERIFIED_DISAGREEMENT") process.exitCode = 3;
}

main(process.argv.slice(2)).catch((error) => {
  process.stderr.write(`error: ${error.message}\n`);
  process.exitCode = 1;
});
