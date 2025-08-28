const INDICATORS = [
	{
		name: "Methyl Orange",
		low: 3.1,
		high: 4.4,
		colors: ["red", "orange", "yellow"],
	},
	{
		name: "Methyl Red",
		low: 4.4,
		high: 6.2,
		colors: ["red", "orange", "yellow"],
	},
	{
		name: "Bromothymol Blue",
		low: 6.0,
		high: 7.6,
		colors: ["yellow", "green", "blue"],
	},
	{
		name: "Phenolphthalein",
		low: 8.2,
		high: 10.0,
		colors: ["colorless", "pink", "fuchsia"],
	},
	{
		name: "Thymol Blue",
		low: 1.2,
		high: 2.8,
		colors: ["red", "orange", "yellow"],
	},
];

const REACTIONS = [
	{
		label: "Strong acid (HCl) vs Strong base (NaOH)",
		analyte: {
			type: "SA",
			C: rand(0.05, 0.12),
			V: rand(20, 35),
		},
		titrant: { type: "SB", C: rand(0.05, 0.12) },
	},
	{
		label: "Weak acid (CH3COOH) vs Strong base (NaOH)",
		analyte: {
			type: "WA",
			Ka: 1.8e-5,
			C: rand(0.05, 0.12),
			V: rand(20, 35),
		},
		titrant: { type: "SB", C: rand(0.05, 0.12) },
	},
	{
		label: "Strong acid (HCl) vs Weak base (NH3)",
		analyte: {
			type: "SA",
			C: rand(0.05, 0.12),
			V: rand(20, 35),
		},
		titrant: { type: "WB", Kb: 1.8e-5, C: rand(0.05, 0.12) },
	},
	{
		label: "Weak acid (HF) vs Weak base (NH3)",
		analyte: {
			type: "WA",
			Ka: 6.8e-4,
			C: rand(0.06, 0.1),
			V: rand(20, 30),
		},
		titrant: { type: "WB", Kb: 1.8e-5, C: rand(0.06, 0.1) },
	},
];

function rand(a, b) {
	return +(a + Math.random() * (b - a)).toFixed(3);
}
function clamp(x, a, b) {
	return Math.max(a, Math.min(b, x));
}

function H(pH) {
	return Math.pow(10, -pH);
}

function pHColor(pH, ind) {
	// Map pH through indicator range to color stops
	if (ind.name === "Phenolphthalein" && pH < 8.2) return "#dbeafe"; // colorless-ish
	const t = (pH - ind.low) / (ind.high - ind.low);
	if (t <= 0) return pickColor(ind.colors[0]);
	if (t >= 1)
		return pickColor(ind.colors[2] || ind.colors[1] || ind.colors[0]);
	// interpolate between colors[0] -> colors[1] -> colors[2]
	const mid = ind.colors[1] || ind.colors[0];
	if (t < 0.5) return mix(pickColor(ind.colors[0]), pickColor(mid), t * 2);
	return mix(pickColor(mid), pickColor(ind.colors[2] || mid), (t - 0.5) * 2);
}
function pickColor(name) {
	const m = {
		red: "#ef4444",
		orange: "#f59e0b",
		yellow: "#facc15",
		green: "#22c55e",
		blue: "#3b82f6",
		pink: "#f472b6",
		fuchsia: "#d946ef",
		colorless: "#dbeafe",
	};
	return m[name] || "#60a5fa";
}
function mix(hex1, hex2, t) {
	const a = parseInt(hex1.slice(1), 16),
		b = parseInt(hex2.slice(1), 16);
	const r = ((a >> 16) & 255) * (1 - t) + ((b >> 16) & 255) * t;
	const g = ((a >> 8) & 255) * (1 - t) + ((b >> 8) & 255) * t;
	const bl = (a & 255) * (1 - t) + (b & 255) * t;
	return `#${(((1 << 24) + (r << 16) + (g << 8) + bl) | 0)
		.toString(16)
		.slice(1)}`;
}

function calcPH(model, Vadd) {
	// V in mL, C in mol/L. Convert to L for moles.
	const Ca = model.analyte.C,
		Va = model.analyte.V / 1000; // L
	const Ct = model.titrant.C,
		Vt = Vadd / 1000; // L
	const nA = Ca * Va; // initial moles analyte (acid or base depending on type)
	const nT = Ct * Vt; // moles titrant added

	// Determine net after neutralization for 4 classes
	// Return pH using approximations
	const typeA = model.analyte.type; // SA or WA
	const typeT = model.titrant.type; // SB or WB

	// Equivalence volume (for reveal after round)
	model.eqV = (nA / Ct) * 1000; // mL

	// Strong-Strong
	if (typeA === "SA" && typeT === "SB") {
		if (nT < nA) {
			// acid in excess
			const Hm = (nA - nT) / (Va + Vt);
			return -Math.log10(Hm);
		} else if (nT > nA) {
			const OHm = (nT - nA) / (Va + Vt);
			return 14 + Math.log10(OHm);
		} else return 7.0;
	}

	// Weak acid + Strong base (buffer/Henderson-Hasselbalch before eq)
	if (typeA === "WA" && typeT === "SB") {
		const Ka = model.analyte.Ka;
		const pKa = -Math.log10(Ka);
		if (nT === 0) {
			// pure weak acid
			const H = Math.sqrt(Ka * Ca);
			return -Math.log10(H);
		}
		if (nT < nA) {
			const ratio = nT / (nA - nT); // base/acid
			return pKa + Math.log10(ratio);
		} else if (nT > nA) {
			const OHm = (nT - nA) / (Va + Vt);
			return 14 + Math.log10(OHm);
		} else {
			// at equivalence: basic from A- hydrolysis
			const Csalt = nA / (Va + Vt);
			const Kb = 1e-14 / Ka;
			const OH = Math.sqrt(Kb * Csalt);
			return 14 + Math.log10(OH);
		}
	}

	// Strong acid + Weak base (mirror)
	if (typeA === "SA" && typeT === "WB") {
		const Kb = model.titrant.Kb;
		const pKb = -Math.log10(Kb);
		const pKa = 14 - pKb;
		if (nT === 0) {
			return 0.5 * (-Math.log10(Ca) - Math.log10(1));
		} // ~very acidic
		if (nT < nA) {
			const ratio = (nA - nT) / nT; // acid/base (use pH = pKa + log([A-]/[HA]) for conjugate acid)
			return pKa - Math.log10(ratio);
		} else if (nT > nA) {
			const Hm = (nT - nA) / (Va + Vt);
			return -Math.log10(Hm);
		} else {
			// equivalence acidic
			const Csalt = nA / (Va + Vt);
			const Ka = 1e-14 / Kb;
			const H = Math.sqrt(Ka * Csalt);
			return -Math.log10(H);
		}
	}

	// Weak acid vs Weak base (approximate using Henderson around midpoint)
	if (typeA === "WA" && typeT === "WB") {
		const Ka = model.analyte.Ka,
			Kb = model.titrant.Kb;
		const pKa = -Math.log10(Ka),
			pKb = -Math.log10(Kb);
		if (nT === 0) {
			const H = Math.sqrt(Ka * Ca);
			return -Math.log10(H);
		}
		if (nT === nA) {
			return 7 + 0.5 * (pKa - pKb);
		}
		if (nT < nA) {
			// buffer region towards acid
			const ratio = nT / (nA - nT);
			return pKa + Math.log10(ratio);
		} else {
			// base excess
			const ratio = (nT - nA) / nA; // crude
			return 14 - (pKb + Math.log10(ratio + 1e-9));
		}
	}

	return 7;
}

// --- UI & Game State ---
const state = {
	round: 1,
	model: null,
	indicator: null,
	players: [
		{
			id: 1,
			name: "Player 1",
			added: 0,
			submitted: false,
			guess: null,
			result: null,
		},
		{
			id: 2,
			name: "Player 2",
			added: 0,
			submitted: false,
			guess: null,
			result: null,
		},
	],
	volumesHidden: false,
};

const grid = document.getElementById("grid");
const newRoundBtn = document.getElementById("newRoundBtn");
const fsBtn = document.getElementById("fsBtn");

fsBtn.onclick = () => {
	const elem = document.documentElement;
	if (!document.fullscreenElement) {
		elem.requestFullscreen?.();
	} else {
		document.exitFullscreen?.();
	}
};

function newRound() {
	// Randomly pick one reaction and one indicator
	const r = JSON.parse(
		JSON.stringify(REACTIONS[Math.floor(Math.random() * REACTIONS.length)])
	);
	// fresh randomization for each round (since REACTIONS were templated with initial rand)
	r.analyte.C = rand(0.05, 0.12);
	r.analyte.V = rand(20, 35);
	r.titrant.C = rand(0.05, 0.12);
	state.model = r;
	state.indicator = INDICATORS[Math.floor(Math.random() * INDICATORS.length)];

	state.players.forEach((p) => {
		p.added = 0;
		p.submitted = false;
		p.guess = null;
		p.result = null;
	});
	state.winner = null;
	state.volumesHidden = false;
	newRoundBtn.disabled = true;
	render();
}

function addVolume(pid, mL) {
	const p = state.players.find((x) => x.id === pid);
	if (p.submitted) return; // can't add after submit
	p.added = +(p.added + mL).toFixed(3);
	render();
}

function submit(pid) {
	const p = state.players.find((x) => x.id === pid);
	if (p.submitted) return;
	p.submitted = true;
	p.guess = p.added;

	// Mark that volumes should be hidden for opponents
	state.volumesHidden = true;

	// Each player can always see their own value
	state.players.forEach((player) => {
		player.showOwn = player.submitted || player.id === pid;
	});

	// Once both submitted, unhide all volumes and decide winner
	if (state.players.every((x) => x.submitted)) {
		state.volumesHidden = false;
		state.players.forEach((player) => (player.showOwn = true));
		decideWinner();
	}

	render();
}

function decideWinner() {
	const eq = state.model.eqV ?? (calcPH(state.model, 0) && state.model.eqV); // ensure eqV set
	state.players.forEach((p) => {
		const err = Math.abs((p.guess ?? 0) - eq);
		p.result = { error: err, eq: eq };
	});
	const [p1, p2] = state.players;
	let winner = null;
	if (p1.result.error === p2.result.error) winner = null;
	else winner = p1.result.error < p2.result.error ? p1.id : p2.id;
	state.winner = winner;
	state.volumesHidden = false; // reveal after both submit
	newRoundBtn.disabled = false;
}

function fmt(x, dec = 2) {
	return (Math.round(x * 10 ** dec) / 10 ** dec).toFixed(dec);
}

function buildPlayerCard(p) {
	const m = state.model,
		ind = state.indicator;
	const pH = clamp(calcPH(m, p.added), 0, 14);
	const color = pHColor(pH, ind);

	const hideVol = state.volumesHidden || false;

	// SVG Flask with colored solution level based on analyte volume + titrant volume (height proportion)
	const totalVol = m.analyte.V + p.added; // mL
	const level = clamp(totalVol / 80, 0.15, 0.95); // visual level (normalized)

	return `
				<div class="card ${
					state.winner
						? state.winner === p.id
							? "winner"
							: "loser"
						: ""
				}">
					<div class="row" style="justify-content:space-between;align-items:flex-start">
					<h2>${p.name}</h2>
					<div class="legend">
						<span class="badge">${m.label}</span>
						<span class="badge">Indicator: ${ind.name} (${ind.low}-${ind.high})</span>
					</div>
					</div>

					<div class="flaskBox">
					<svg class="flask" viewBox="0 0 120 180" aria-label="conical flask">
						<!-- Flask outline -->
						<defs>
						<linearGradient id="g${p.id}" x1="0" x2="0" y1="0" y2="1">
							<stop offset="0%" stop-color="${color}" stop-opacity="0.9"/>
							<stop offset="100%" stop-color="${color}" stop-opacity="0.95"/>
						</linearGradient>
						</defs>
						<path d="M45 10 h30 v60 c0 8 3 10 10 18 l20 40 c5 10 -2 20 -12 20 H32 c-10 0 -17 -10 -12 -20 l20 -40 c7 -8 10 -10 10 -18 V10z" fill="#0b1220" stroke="#334155" stroke-width="2"/>
						<!-- Liquid -->
						<clipPath id="clip${p.id}">
						<path d="M45 10 h30 v60 c0 8 3 10 10 18 l20 40 c5 10 -2 20 -12 20 H32 c-10 0 -17 -10 -12 -20 l20 -40 c7 -8 10 -10 10 -18 V10z"/>
						</clipPath>
						<g clip-path="url(#clip${p.id})">
						<rect x="10" y="${180 - 150 * level}" width="100" height="150" fill="url(#g${
		p.id
	})"/>
						<!-- tiny bubbles -->
						<circle cx="50" cy="${170 - 150 * level}" r="1.5" fill="#ffffff22"/>
						<circle cx="65" cy="${160 - 150 * level}" r="1" fill="#ffffff22"/>
						<circle cx="80" cy="${150 - 150 * level}" r="1.2" fill="#ffffff22"/>
						</g>
					</svg>
					<div>
						<div class="row">
						<div class="stat">
							<label>pH (simulated)</label>
							<div class="val">${fmt(pH, 2)}</div>
						</div>
						<div class="stat">
							<label>Titrant added (mL)</label>
							<div class="val ${
								state.volumesHidden && p.submitted
									? "blurred"
									: ""
							}" title="Hidden until both submit">${fmt(
		p.added,
		2
	)}</div>
						</div>
						<div class="stat">
							<label>Analyte Vol (mL)</label>
							<div class="val">${fmt(m.analyte.V, 1)}</div>
						</div>
						<div class="stat">
							<label>Conc (analyte/titrant, M)</label>
							<div class="val">${fmt(m.analyte.C, 3)} / ${fmt(m.titrant.C, 3)}</div>
						</div>
						</div>
						<div class="meter" aria-label="progress to equivalence (hidden)"><div style="width:${clamp(
							(p.added / (state.model.eqV || 100)) * 100,
							0,
							100
						)}%"></div></div>
						<div class="controlsBar" style="margin-top:10px">
						<button onclick="addVolume(${p.id},0.1)" ${
		p.submitted ? "disabled" : ""
	}>+0.10 mL</button>
						<button onclick="addVolume(${p.id},0.5)" ${
		p.submitted ? "disabled" : ""
	}>+0.50 mL</button>
						<button onclick="addVolume(${p.id},1.0)" ${
		p.submitted ? "disabled" : ""
	}>+1.00 mL</button>
						<button onclick="addVolume(${p.id},2.0)" ${
		p.submitted ? "disabled" : ""
	}>+2.00 mL</button>
						<button onclick="submit(${p.id})" ${
		p.submitted ? "disabled" : ""
	}>Submit</button>
						</div>
						${
							p.submitted && p.result
								? `<div class="row" style="margin-top:8px;gap:14px">
						<span class="pill">Your guess: ${fmt(p.guess, 2)} mL</span>
						<span class="pill">Equivalence: ${fmt(p.result.eq, 2)} mL</span>
						<span class="pill">Absolute error: ${fmt(p.result.error, 2)} mL</span>
						</div>`
								: ""
						}
					</div>
					</div>
				</div>`;
}

function render() {
	// Ensure equivalence cached
	calcPH(state.model, 0);
	grid.innerHTML = state.players.map(buildPlayerCard).join("");

	// If both submitted, show banner
	const both = state.players.every((p) => p.submitted);
	const notice = document.getElementById("notice");
	if (both) {
		if (state.winner) {
			notice.innerHTML = `<b>Round ${state.round}:</b> Winner is Player ${state.winner}! Start next round when ready.`;
		} else {
			notice.innerHTML = `<b>Round ${state.round}:</b> It's a tie! Start next round when ready.`;
		}
	} else {
		const who = state.players
			.filter((p) => p.submitted)
			.map((p) => p.id)
			.join(", ");
		notice.innerHTML = who
			? `Player ${who} submitted. Volumes hidden until both submit.`
			: `By Justin and Arav`;
	}
}

// Keyboard shortcuts for quick play
window.addEventListener("keydown", (e) => {
	const k = e.key.toLowerCase();
	if (k === "a") addVolume(1, 0.1);
	if (k === "s") addVolume(1, 0.5);
	if (k === "d") addVolume(1, 1.0);
	if (k === "f") submit(1);

	if (k === "j") addVolume(2, 0.1);
	if (k === "k") addVolume(2, 0.5);
	if (k === "l") addVolume(2, 1.0);
	if (k === ";") submit(2);
});

newRoundBtn.addEventListener("click", () => {
	state.round++;
	newRound();
});

// Init
newRound();
