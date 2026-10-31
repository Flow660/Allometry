# app.py
# Streamlitアプリ（全コード）
# ------------------------------------------------------------
# 目的（ユーザー要件）:
#   1) 3つの予測パターン（2種→残り1種）で
#        - Mean（2種平均）
#        - 2 species allometry（2点回帰: t=a*BW^b）
#      を比較し、各パターンで「どちらがベストか」を出力する。
#
#      パターン:
#        A) Mouse + Rat    -> Monkey
#        B) Mouse + Monkey -> Rat
#        C) Rat + Monkey   -> Mouse
#
#   2) 追加要件:
#      ① 2種allometryで求めた定数（a,b）を同じCSVに出力する
#      ② 「3種allometry（3点回帰）」と「3種平均（全3種平均）」で
#         各動物種（Mouse/Rat/Monkey）を予測した場合に、どちらが精度が高いかの比較を追加。
#         ②の列は、これまでのCSVの右側にまとめて追加する。
#
# 入力:
#   Excel (.xlsx) に 4列:
#     Protein名 / Mouse半減期 / Rat半減期 / Monkey半減期
#   ※列名が多少違っても、先頭4列を使うフォールバックあり
#
# 体重(kg)（固定）:
#   Mouse=0.025, Rat=0.2, Monkey=6
#
# 誤差（予測精度指標）:
#   倍率誤差を公平に扱うため log空間の誤差を使用:
#     - 個別:  ErrorAbs = |ln(true) - ln(pred)|
#     - 全体:  ErrorRMSE = sqrt(mean( (ln(true)-ln(pred))^2 ))
#
# 注意:
#   添付Excelには欠損や0/負値が含まれるため、0以下/欠損はlogが取れず比較不能。
#   その場合は該当の予測/誤差/Best判定を NA（比較不可）にします。
#
# 実行:
#   pip install streamlit pandas numpy openpyxl matplotlib   (scipy は任意: あれば Wilcoxon 検定に使用)
#   streamlit run app.py
#
# 追加（Figure/Table 出力）:
#   4) Figure: (a) Mean 予測 vs 観測, (b) 2点Allometry 予測 vs 観測, (c) 箱ひげ図 + Wilcoxon検定
#      Table 1: 3パターンで Mean / 2-species allometry が優れたタンパク質数
#      いずれも PNG でダウンロード可能
#
# 追加（比較用の基準 / Baseline, 2026-10-27 v3）:
#   5) 予測ベンチマークの比較基準（ターゲット種のデータは一切使わない）
#      ① 一定値の予測: 全タンパク質に同じ値 = 2種平均予測の全タンパク質での中央値
#      ② 指数固定の allometry: t_target = 幾何平均( t_s1*(BW_T/BW_s1)^b , t_s2*(BW_T/BW_s2)^b ), b = 0.25（既定）
#         （生理学的時間のスケーリング: Boxenbaum, J Pharmacokinet Biopharm 1982）
#      指標: 2倍/3倍以内(%), GMFE, RMSLE(log10), Pearson r(log10), |log10 誤差|の中央値
#      CSV には各パターンの予測値・誤差列（Pred*Const*, Pred*FixedAllo*）を追加
# ------------------------------------------------------------

import io
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import streamlit as st


# =========================
# 固定パラメータ（要件）
# =========================
BW = {
    "Mouse": 0.025,
    "Rat": 0.2,
    "Monkey": 6.0,
}

DEFAULT_EXCEL = "ProteinList_common.xlsx"
DEFAULT_SHEET = "CommonList"   # 3種すべての値が揃ったシート


# =========================
# 基本ユーティリティ
# =========================
def is_valid_pos(x) -> bool:
    """正の実数ならTrue（logが取れる）"""
    try:
        return (x is not None) and (not pd.isna(x)) and (float(x) > 0.0)
    except Exception:
        return False


def ln(x: float) -> float:
    return float(np.log(float(x)))


def pred_mean(x1: float, x2: float) -> float:
    return (float(x1) + float(x2)) / 2.0


def fit_allometry_two_points(t1: float, w1: float, t2: float, w2: float):
    """
    2点で allometry を作る:
      ln(t)=ln(a)+b ln(BW)
    戻り値: (a, b)
    """
    lnw1, lnw2 = np.log(w1), np.log(w2)
    lnt1, lnt2 = np.log(t1), np.log(t2)

    b = (lnt2 - lnt1) / (lnw2 - lnw1)
    lna = lnt1 - b * lnw1
    a = float(np.exp(lna))
    return a, float(b)


def pred_allometry(a: float, b: float, wT: float) -> float:
    """t = a * BW^b"""
    return float(a) * (float(wT) ** float(b))


def fit_allometry_three_points(tM: float, tR: float, tK: float):
    """
    3点（Mouse/Rat/Monkey）で OLS 回帰し allometry をフィット:
      y = ln(t), x = ln(BW)
      y = alpha + b x  (alpha=ln(a))
    戻り値: (a, b)
    """
    xs = np.array([np.log(BW["Mouse"]), np.log(BW["Rat"]), np.log(BW["Monkey"])], dtype=float)
    ys = np.array([np.log(tM), np.log(tR), np.log(tK)], dtype=float)

    xbar = xs.mean()
    ybar = ys.mean()
    denom = np.sum((xs - xbar) ** 2)
    if denom == 0:
        return np.nan, np.nan

    b = float(np.sum((xs - xbar) * (ys - ybar)) / denom)
    alpha = float(ybar - b * xbar)  # = ln(a)
    a = float(np.exp(alpha))
    return a, b


def err_log_abs(true_t: float, pred_t: float) -> float:
    """|ln(true)-ln(pred)|（倍率誤差）"""
    if (not is_valid_pos(true_t)) or (not is_valid_pos(pred_t)):
        return np.nan
    return float(abs(ln(true_t) - ln(pred_t)))


def err_log_rmse(true_list, pred_list) -> float:
    """RMSE in log-space: sqrt(mean((ln(true)-ln(pred))^2))"""
    pairs = []
    for t, p in zip(true_list, pred_list):
        if is_valid_pos(t) and is_valid_pos(p):
            pairs.append((ln(t) - ln(p)) ** 2)
    if not pairs:
        return np.nan
    return float(np.sqrt(np.mean(pairs)))


def best_of_two(err_mean: float, err_allo: float) -> str:
    """2方式のうち誤差が小さい方を返す（同点は同等）"""
    if pd.isna(err_mean) and pd.isna(err_allo):
        return "比較不可(NA)"
    if pd.isna(err_mean):
        return "Allometry"
    if pd.isna(err_allo):
        return "Mean"

    if err_mean < err_allo:
        return "Mean"
    if err_allo < err_mean:
        return "Allometry"
    return "同等(誤差同じ)"


# =========================
# 入力列の正規化
# =========================
def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    入力DFの列名をそれっぽく解釈して統一する。
    期待する出力列: Protein, Mouse, Rat, Monkey
    """
    cols = list(df.columns)
    lower = {c: str(c).strip().lower() for c in cols}

    def find_col(candidates):
        for c in cols:
            if lower[c] in candidates:
                return c
        return None

    protein_col = find_col({"protein", "proteinname", "pg.proteingroups", "protein groups", "proteingroups", "name"})
    mouse_col = find_col({"mouse", "mice"})
    rat_col = find_col({"rat"})
    monkey_col = find_col({"monkey", "monkeys", "primate"})

    # 見つからない場合は「先頭4列」を採用（A-D形式に対応）
    if protein_col is None or mouse_col is None or rat_col is None or monkey_col is None:
        if len(cols) >= 4:
            protein_col, mouse_col, rat_col, monkey_col = cols[0], cols[1], cols[2], cols[3]
        else:
            raise ValueError("入力に必要な4列(Protein/Mouse/Rat/Monkey)が見つかりませんでした。")

    out = df.rename(
        columns={
            protein_col: "Protein",
            mouse_col: "Mouse",
            rat_col: "Rat",
            monkey_col: "Monkey",
        }
    ).copy()

    for c in ["Mouse", "Rat", "Monkey"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    out["Protein"] = out["Protein"].astype(str)
    return out[["Protein", "Mouse", "Rat", "Monkey"]]


# =========================
# 1行（1Protein）計算
# =========================
def compute_row(row: pd.Series) -> dict:
    """
    左側: 3パターン（2種→1種）で Mean vs 2点Allometry を比較 + 2点定数(a,b)出力
    右側: 3種Mean vs 3点Allometry（3点回帰）比較（Mouse/Rat/Monkeyそれぞれ + 全体）
    """
    tM, tR, tK = row["Mouse"], row["Rat"], row["Monkey"]  # MonkeyをKと表記

    out_left = {}

    # =========================================================
    # A) Mouse+Rat -> Monkey
    # =========================================================
    # Allometry2定数は「source2点が正」なら出す（targetが無くてもOK）
    if is_valid_pos(tM) and is_valid_pos(tR):
        a_MR, b_MR = fit_allometry_two_points(tM, BW["Mouse"], tR, BW["Rat"])
        predA_allo = pred_allometry(a_MR, b_MR, BW["Monkey"])
        predA_mean = pred_mean(tM, tR)
    else:
        a_MR = b_MR = np.nan
        predA_allo = predA_mean = np.nan

    errA_mean = err_log_abs(tK, predA_mean)
    errA_allo = err_log_abs(tK, predA_allo)
    bestA = best_of_two(errA_mean, errA_allo)

    out_left.update({
        # constants
        "Allo2_a_MouseRat": a_MR,
        "Allo2_b_MouseRat": b_MR,
        # predictions / errors / best
        "PredMonkeyMeanFromMouseRat": predA_mean,
        "PredMonkeyAllo2FromMouseRat": predA_allo,
        "ErrMonkeyMeanFromMouseRat": errA_mean,
        "ErrMonkeyAllo2FromMouseRat": errA_allo,
        "BestMonkey_MeanVsAllo2": bestA,
    })

    # =========================================================
    # B) Mouse+Monkey -> Rat
    # =========================================================
    if is_valid_pos(tM) and is_valid_pos(tK):
        a_MK, b_MK = fit_allometry_two_points(tM, BW["Mouse"], tK, BW["Monkey"])
        predB_allo = pred_allometry(a_MK, b_MK, BW["Rat"])
        predB_mean = pred_mean(tM, tK)
    else:
        a_MK = b_MK = np.nan
        predB_allo = predB_mean = np.nan

    errB_mean = err_log_abs(tR, predB_mean)
    errB_allo = err_log_abs(tR, predB_allo)
    bestB = best_of_two(errB_mean, errB_allo)

    out_left.update({
        "Allo2_a_MouseMonkey": a_MK,
        "Allo2_b_MouseMonkey": b_MK,
        "PredRatMeanFromMouseMonkey": predB_mean,
        "PredRatAllo2FromMouseMonkey": predB_allo,
        "ErrRatMeanFromMouseMonkey": errB_mean,
        "ErrRatAllo2FromMouseMonkey": errB_allo,
        "BestRat_MeanVsAllo2": bestB,
    })

    # =========================================================
    # C) Rat+Monkey -> Mouse
    # =========================================================
    if is_valid_pos(tR) and is_valid_pos(tK):
        a_RK, b_RK = fit_allometry_two_points(tR, BW["Rat"], tK, BW["Monkey"])
        predC_allo = pred_allometry(a_RK, b_RK, BW["Mouse"])
        predC_mean = pred_mean(tR, tK)
    else:
        a_RK = b_RK = np.nan
        predC_allo = predC_mean = np.nan

    errC_mean = err_log_abs(tM, predC_mean)
    errC_allo = err_log_abs(tM, predC_allo)
    bestC = best_of_two(errC_mean, errC_allo)

    out_left.update({
        "Allo2_a_RatMonkey": a_RK,
        "Allo2_b_RatMonkey": b_RK,
        "PredMouseMeanFromRatMonkey": predC_mean,
        "PredMouseAllo2FromRatMonkey": predC_allo,
        "ErrMouseMeanFromRatMonkey": errC_mean,
        "ErrMouseAllo2FromRatMonkey": errC_allo,
        "BestMouse_MeanVsAllo2": bestC,
        "All3Positive": bool(is_valid_pos(tM) and is_valid_pos(tR) and is_valid_pos(tK)),
    })

    # =========================================================
    # 右側: ② 3種Mean vs 3点Allometry（3点回帰）比較
    # =========================================================
    out_right = {}

    if is_valid_pos(tM) and is_valid_pos(tR) and is_valid_pos(tK):
        # Mean3
        mean3 = float((tM + tR + tK) / 3.0)
        predM_mean3 = mean3
        predR_mean3 = mean3
        predK_mean3 = mean3

        # Allometry3 fit
        a3, b3 = fit_allometry_three_points(tM, tR, tK)
        predM_allo3 = pred_allometry(a3, b3, BW["Mouse"])
        predR_allo3 = pred_allometry(a3, b3, BW["Rat"])
        predK_allo3 = pred_allometry(a3, b3, BW["Monkey"])

        # Per-species errors (abs log error)
        eM_mean3 = err_log_abs(tM, predM_mean3)
        eR_mean3 = err_log_abs(tR, predR_mean3)
        eK_mean3 = err_log_abs(tK, predK_mean3)

        eM_allo3 = err_log_abs(tM, predM_allo3)
        eR_allo3 = err_log_abs(tR, predR_allo3)
        eK_allo3 = err_log_abs(tK, predK_allo3)

        bestM = best_of_two(eM_mean3, eM_allo3)
        bestR = best_of_two(eR_mean3, eR_allo3)
        bestK = best_of_two(eK_mean3, eK_allo3)

        # Overall RMSE in log space
        rmse_mean3 = err_log_rmse([tM, tR, tK], [predM_mean3, predR_mean3, predK_mean3])
        rmse_allo3 = err_log_rmse([tM, tR, tK], [predM_allo3, predR_allo3, predK_allo3])
        best_overall = best_of_two(rmse_mean3, rmse_allo3)

    else:
        mean3 = np.nan
        a3 = b3 = np.nan
        predM_mean3 = predR_mean3 = predK_mean3 = np.nan
        predM_allo3 = predR_allo3 = predK_allo3 = np.nan
        eM_mean3 = eR_mean3 = eK_mean3 = np.nan
        eM_allo3 = eR_allo3 = eK_allo3 = np.nan
        bestM = bestR = bestK = "比較不可(NA)"
        rmse_mean3 = rmse_allo3 = np.nan
        best_overall = "比較不可(NA)"

    out_right.update({
        # 3species mean
        "Mean3_AllSpecies": mean3,

        # 3species allometry constants
        "Allo3_a_AllSpecies": a3,
        "Allo3_b_AllSpecies": b3,

        # predictions from mean3 / allo3
        "PredMouse_Mean3": predM_mean3,
        "PredMouse_Allo3": predM_allo3,
        "PredRat_Mean3": predR_mean3,
        "PredRat_Allo3": predR_allo3,
        "PredMonkey_Mean3": predK_mean3,
        "PredMonkey_Allo3": predK_allo3,

        # errors per species
        "ErrMouse_Mean3": eM_mean3,
        "ErrMouse_Allo3": eM_allo3,
        "BestMouse_Mean3VsAllo3": bestM,

        "ErrRat_Mean3": eR_mean3,
        "ErrRat_Allo3": eR_allo3,
        "BestRat_Mean3VsAllo3": bestR,

        "ErrMonkey_Mean3": eK_mean3,
        "ErrMonkey_Allo3": eK_allo3,
        "BestMonkey_Mean3VsAllo3": bestK,

        # overall comparison (RMSE log)
        "RMSElog_Mean3_AllSpecies": rmse_mean3,
        "RMSElog_Allo3_AllSpecies": rmse_allo3,
        "BestOverall_Mean3VsAllo3": best_overall,
    })

    # 順序を「左(既存)→右(追加)」にする
    out = {}
    out.update(out_left)
    out.update(out_right)
    return out


def run_analysis(df_in: pd.DataFrame) -> pd.DataFrame:
    df = normalize_columns(df_in)
    extra = df.apply(lambda r: compute_row(r), axis=1, result_type="expand")
    return pd.concat([df, extra], axis=1)


# =========================
# 論文用 Figure / Table 出力（追加）
# =========================
# パターン定義: (source1, source2, target, 表示ラベル, Mean予測列, Allo2予測列, Mean誤差列, Allo2誤差列)
PATTERNS = {
    "Monkey": ("Mouse", "Rat", "Monkey", "Mouse&Rat\u2192Monkey",
               "PredMonkeyMeanFromMouseRat", "PredMonkeyAllo2FromMouseRat",
               "ErrMonkeyMeanFromMouseRat", "ErrMonkeyAllo2FromMouseRat"),
    "Rat": ("Mouse", "Monkey", "Rat", "Mouse&Monkey\u2192Rat",
            "PredRatMeanFromMouseMonkey", "PredRatAllo2FromMouseMonkey",
            "ErrRatMeanFromMouseMonkey", "ErrRatAllo2FromMouseMonkey"),
    "Mouse": ("Rat", "Monkey", "Mouse", "Rat&Monkey\u2192Mouse",
              "PredMouseMeanFromRatMonkey", "PredMouseAllo2FromRatMonkey",
              "ErrMouseMeanFromRatMonkey", "ErrMouseAllo2FromRatMonkey"),
}

FIXED_EXPONENT_DEFAULT = 0.25   # 生理学的時間のスケーリング（t ∝ BW^0.25）


def _pattern_keys(target: str):
    s1, s2, t = PATTERNS[target][:3]
    return s1, s2, t, f"{s1}{s2}"


def add_baselines(df_out: pd.DataFrame, exponent: float = FIXED_EXPONENT_DEFAULT) -> pd.DataFrame:
    """
    比較用の基準（baseline）の予測値と誤差 |ln(true)-ln(pred)| を各パターンについて追加する。
      ① 一定値: 2種平均予測（source 2種の算術平均）の、解析対象全タンパク質での中央値を全タンパク質に適用
      ② 指数固定 allometry: source 2種それぞれを (BW_T/BW_s)^b で補正し、その幾何平均を予測値とする
    いずれもターゲット種の観測値は使わない。
    """
    df = df_out.copy()
    for target in ["Monkey", "Rat", "Mouse"]:
        s1, s2, t, tag = _pattern_keys(target)
        c_pm = PATTERNS[target][4]
        pm = pd.to_numeric(df[c_pm], errors="coerce")
        const = float(np.nanmedian(pm[pm > 0])) if (pm > 0).any() else np.nan
        pc = np.where(pm > 0, const, np.nan)
        v1 = pd.to_numeric(df[s1], errors="coerce"); v2 = pd.to_numeric(df[s2], errors="coerce")
        ok = (v1 > 0) & (v2 > 0)
        f1 = (BW[t] / BW[s1]) ** exponent; f2 = (BW[t] / BW[s2]) ** exponent
        pf = np.where(ok, np.sqrt(v1.where(ok, 1.0) * f1 * v2.where(ok, 1.0) * f2), np.nan)
        df[f"Pred{t}ConstFrom{tag}"] = pc
        df[f"Pred{t}FixedAlloFrom{tag}"] = pf
        df[f"Err{t}ConstFrom{tag}"] = [err_log_abs(o, p) for o, p in zip(df[t], pc)]
        df[f"Err{t}FixedAlloFrom{tag}"] = [err_log_abs(o, p) for o, p in zip(df[t], pf)]
        df.attrs[f"Const_{t}"] = const
    df.attrs["FixedExponent"] = exponent
    return df


def _metrics(obs, pred) -> dict:
    lo, lp = np.log10(np.asarray(obs, float)), np.log10(np.asarray(pred, float))
    e = lp - lo; a = np.abs(e)
    r = float(np.corrcoef(lp, lo)[0, 1]) if np.ptp(lp) > 1e-9 else np.nan   # 一定値の予測では未定義
    return {
        "Within 2-fold (%)": float((a <= np.log10(2)).mean() * 100),
        "Within 3-fold (%)": float((a <= np.log10(3)).mean() * 100),
        "GMFE": float(10 ** a.mean()),
        "RMSLE (log10)": float(np.sqrt((e ** 2).mean())),
        "Pearson r (log10)": r,
        "Median |log10 error|": float(np.median(a)),
    }


def make_baseline_table(df_b: pd.DataFrame, target: str = "Monkey") -> pd.DataFrame:
    """
    補足表: 同じタンパク質集合（全方式の予測が得られた行）で 4方式を比較
      Mean / 2-species allometry / ① 一定値 / ② 指数固定 allometry
    """
    s1, s2, t, tag = _pattern_keys(target)
    c_pm, c_pa = PATTERNS[target][4], PATTERNS[target][5]
    cols = {"Mean (species mean)": c_pm, "2-species allometry": c_pa,
            "Constant (median of species-mean predictions)": f"Pred{t}ConstFrom{tag}",
            f"Fixed-exponent allometry (b = {df_b.attrs.get('FixedExponent', FIXED_EXPONENT_DEFAULT):g})": f"Pred{t}FixedAlloFrom{tag}"}
    obs = pd.to_numeric(df_b[t], errors="coerce")
    ok = obs > 0
    for c in cols.values():
        ok &= pd.to_numeric(df_b[c], errors="coerce") > 0
    rows = {}
    for lab, c in cols.items():
        m = _metrics(obs[ok], pd.to_numeric(df_b.loc[ok, c]))
        rows[lab] = {"n": int(ok.sum()), **m}
    tbl = pd.DataFrame(rows).T
    tbl.index.name = PATTERNS[target][3]
    return tbl


def make_baseline_table_figure(tbl: pd.DataFrame, const_value: float, exponent: float, dpi: int = 300):
    """補足表（比較用の基準）を論文風の PNG 用 Figure として描画"""
    cols = ["n", "Within 2-fold (%)", "Within 3-fold (%)", "GMFE", "RMSLE (log10)", "Pearson r (log10)"]
    fmt = {"n": "{:,.0f}", "Within 2-fold (%)": "{:.1f}", "Within 3-fold (%)": "{:.1f}", "GMFE": "{:.2f}",
           "RMSLE (log10)": "{:.2f}", "Pearson r (log10)": "{:.2f}"}
    row_in = 0.42; fig_w = 12.5; fig_h = 2.2 + row_in * (len(tbl) + 1)
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, fig_w); ax.set_ylim(0, fig_h)
    y = fig_h - 0.35
    ax.text(0.2, y, f"Prediction of {tbl.index.name}: species mean and allometry versus reference predictors",
            fontsize=11, fontweight="bold", va="center")
    x_lab = 0.25; col_x = [5.3, 6.4, 7.6, 8.7, 9.8, 11.2]
    x0, x1 = 0.2, fig_w - 0.2
    y -= 0.55; ax.plot([x0, x1], [y, y], color="black", lw=1.3)
    y -= row_in / 2
    ax.text(x_lab, y, "Method", ha="left", va="center", fontsize=10)
    for x, c in zip(col_x, cols):
        ax.text(x, y, c, ha="center", va="center", fontsize=9.5)
    y -= row_in / 2; ax.plot([x0, x1], [y, y], color="black", lw=0.8)
    for r in tbl.index:
        y -= row_in / 2
        ax.text(x_lab, y, r, ha="left", va="center", fontsize=10)
        for x, c in zip(col_x, cols):
            v = tbl.loc[r, c]
            ax.text(x, y, "—" if pd.isna(v) else fmt[c].format(v), ha="center", va="center", fontsize=10)
        y -= row_in / 2
    ax.plot([x0, x1], [y, y], color="black", lw=1.3)
    y -= 0.3
    ax.text(0.2, y,
            f"Constant: every protein assigned the median species-mean prediction ({const_value:.2f} days). "
            f"Fixed-exponent allometry: geometric mean of source half-lives scaled by (BW_target/BW_source)^{exponent:g}.\n"
            "Neither reference uses target-species data. GMFE, geometric mean fold-error; RMSLE, root-mean-square log10 error; "
            "r is undefined for the constant predictor.",
            fontsize=8.5, va="top")
    return fig


COLOR_MEAN = "#1f77b4"   # 青（Mean）
COLOR_ALLO = "#d62728"   # 赤（Allometry）


def wilcoxon_signed_rank(x, y):
    """
    対応のある2群の両側 Wilcoxon signed-rank 検定。
    scipy があれば scipy.stats.wilcoxon を使用し、無ければ正規近似（同順位補正あり、
    差=0 は除外 [zero_method='wilcox']）で計算する。
    戻り値: (統計量W, p値, 有効ペア数n)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    d = x[ok] - y[ok]
    d = d[d != 0]
    n = int(d.size)
    if n == 0:
        return np.nan, np.nan, 0
    try:
        from scipy.stats import wilcoxon  # type: ignore
        res = wilcoxon(x[ok], y[ok], zero_method="wilcox", alternative="two-sided")
        return float(res.statistic), float(res.pvalue), n
    except Exception:
        pass
    ranks = pd.Series(np.abs(d)).rank(method="average").to_numpy()
    w_plus = float(ranks[d > 0].sum())
    w_minus = float(ranks[d < 0].sum())
    mu = n * (n + 1) / 4.0
    _, counts = np.unique(np.abs(d), return_counts=True)
    tie = float(np.sum(counts ** 3 - counts)) / 48.0
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0 - tie)
    if sigma == 0:
        return min(w_plus, w_minus), np.nan, n
    z = (w_plus - mu) / sigma
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return min(w_plus, w_minus), float(p), n


def format_p(p: float) -> str:
    """P値を論文表記（P < 1 × 10^-k など）に整形"""
    if p is None or pd.isna(p):
        return "P = NA"
    if p == 0.0:
        return "P < 1 \u00d7 10$^{-300}$"
    if p < 1e-3:
        k = int(math.floor(math.log10(p)))
        return f"P < 1 \u00d7 10$^{{{k + 1}}}$"
    return f"P = {p:.3g}"


def pattern_frame(df_out: pd.DataFrame, target: str) -> pd.DataFrame:
    """指定パターンについて、観測値と Mean/Allometry 予測が揃った（正の値の）行のみを返す"""
    s1, s2, t, label, c_pm, c_pa, c_em, c_ea = PATTERNS[target]
    d = pd.DataFrame({
        "Protein": df_out["Protein"],
        "Observed": df_out[t],
        "PredMean": df_out[c_pm],
        "PredAllo": df_out[c_pa],
    })
    ok = d["Observed"].map(is_valid_pos) & d["PredMean"].map(is_valid_pos) & d["PredAllo"].map(is_valid_pos)
    d = d[ok].copy()
    d["Log10RatioMean"] = np.log10(d["PredMean"] / d["Observed"])
    d["Log10RatioAllo"] = np.log10(d["PredAllo"] / d["Observed"])
    d["AbsLog10ErrMean"] = d["Log10RatioMean"].abs()
    d["AbsLog10ErrAllo"] = d["Log10RatioAllo"].abs()
    return d


def _scatter_panel(ax, obs, pred, color, title, lim, show_ylabel):
    """観測 vs 予測の散布図（log-log, 1:1線, 2倍/3倍バンド, 範囲外は三角で端に表示）"""
    lo, hi = lim
    obs = np.asarray(obs, dtype=float)
    pred = np.asarray(pred, dtype=float)

    xx = np.array([lo, hi])
    ax.plot(xx, xx, color="black", lw=1.5, zorder=3)
    for f, ls in [(2.0, "--"), (3.0, ":")]:
        ax.plot(xx, xx * f, color="black", lw=0.9, ls=ls, zorder=3)
        ax.plot(xx, xx / f, color="black", lw=0.9, ls=ls, zorder=3)

    in_x = (obs >= lo) & (obs <= hi)
    inside = in_x & (pred >= lo) & (pred <= hi)
    above = in_x & (pred > hi)
    below = in_x & (pred < lo)

    ax.scatter(obs[inside], pred[inside], s=7, color=color, alpha=0.35, edgecolors="none", zorder=2)
    # 範囲外の予測値は上端/下端に三角で表示
    y_top = hi / 1.06
    y_bot = lo * 1.06
    ax.scatter(obs[above], np.full(above.sum(), y_top), marker="^", s=22,
               facecolors="none", edgecolors=color, linewidths=0.9, zorder=4, clip_on=False)
    ax.scatter(obs[below], np.full(below.sum(), y_bot), marker="^", s=22,
               facecolors="none", edgecolors=color, linewidths=0.9, zorder=4, clip_on=False)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, which="major", color="0.92", lw=0.6)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("Observed half-life (days)", fontsize=9)
    if show_ylabel:
        ax.set_ylabel("Predicted half-life (days)", fontsize=9)
    ax.tick_params(labelsize=8)
    return int(above.sum()), int(below.sum())


def _box_panel(ax, v_mean, v_allo, signed: bool, whisker: float = 1.5):
    """Mean vs Allometry の箱ひげ図（外れ値は非表示、ひげは 1.5×IQR 内の最端値）"""
    bp = ax.boxplot(
        [v_mean, v_allo],
        positions=[1, 2],
        widths=0.55,
        whis=whisker,
        showfliers=False,
        patch_artist=True,
        medianprops=dict(color="black", lw=2.0),
        whiskerprops=dict(color="black", lw=1.0),
        capprops=dict(color="black", lw=1.0),
        boxprops=dict(lw=0.8),
    )
    for patch, col in zip(bp["boxes"], [COLOR_MEAN, COLOR_ALLO]):
        patch.set_facecolor(col)
        patch.set_alpha(0.35)
        patch.set_edgecolor(col)

    l2, l3 = np.log10(2.0), np.log10(3.0)
    ax.axhline(l2, color="0.4", ls="--", lw=0.8)
    ax.axhline(l3, color="0.4", ls=":", lw=0.8)
    if signed:
        ax.axhline(0.0, color="black", lw=1.0)
        ax.axhline(-l2, color="0.4", ls="--", lw=0.8)
        ax.axhline(-l3, color="0.4", ls=":", lw=0.8)
        ax.set_ylabel("log10(Predicted / Observed)", fontsize=9)
    else:
        ax.set_ylabel("|log10(Predicted / Observed)|", fontsize=9)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Mean", "Allometry"], fontsize=9)
    ax.set_xlim(0.4, 2.6)
    ax.tick_params(labelsize=8)


def make_allometry_figure(df_out: pd.DataFrame, target: str = "Monkey",
                          lim=(0.5, 200.0), signed_box: bool = True, dpi: int = 300):
    """
    Figure 5 形式の図を作成:
      (a) Mean 予測 vs 観測, (b) Allometry 予測 vs 観測（同一軸範囲）, (c) 箱ひげ図
    戻り値: (fig, stats dict)
    """
    d = pattern_frame(df_out, target)
    n = len(d)
    _, p, _ = wilcoxon_signed_rank(d["AbsLog10ErrMean"], d["AbsLog10ErrAllo"])

    fig = plt.figure(figsize=(9.0, 9.6), dpi=dpi)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.28, wspace=0.18)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])

    out_a = _scatter_panel(ax_a, d["Observed"], d["PredMean"], COLOR_MEAN, "Mean", lim, True)
    out_b = _scatter_panel(ax_b, d["Observed"], d["PredAllo"], COLOR_ALLO, "Allometry", lim, False)

    if signed_box:
        _box_panel(ax_c, d["Log10RatioMean"], d["Log10RatioAllo"], signed=True)
    else:
        _box_panel(ax_c, d["AbsLog10ErrMean"], d["AbsLog10ErrAllo"], signed=False)
    ax_c.set_title(f"{format_p(p)}  (two-sided Wilcoxon signed-rank on |log10 fold-error|, n = {n:,})",
                   fontsize=8, loc="left")

    for ax, lab in [(ax_a, "a"), (ax_b, "b"), (ax_c, "c")]:
        ax.text(-0.16, 1.04, lab, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")

    stats = {
        "Pattern": PATTERNS[target][3],
        "n": n,
        "Wilcoxon_p_absLog10Err": p,
        "Median_absLog10Err_Mean": float(d["AbsLog10ErrMean"].median()) if n else np.nan,
        "Median_absLog10Err_Allometry": float(d["AbsLog10ErrAllo"].median()) if n else np.nan,
        "Within2fold_Mean_%": float((d["AbsLog10ErrMean"] <= np.log10(2)).mean() * 100) if n else np.nan,
        "Within2fold_Allometry_%": float((d["AbsLog10ErrAllo"] <= np.log10(2)).mean() * 100) if n else np.nan,
        "Within3fold_Mean_%": float((d["AbsLog10ErrMean"] <= np.log10(3)).mean() * 100) if n else np.nan,
        "Within3fold_Allometry_%": float((d["AbsLog10ErrAllo"] <= np.log10(3)).mean() * 100) if n else np.nan,
        "OutOfRange_Mean(above/below)": f"{out_a[0]}/{out_a[1]}",
        "OutOfRange_Allometry(above/below)": f"{out_b[0]}/{out_b[1]}",
    }
    return fig, stats


def make_comparison_table(df_out: pd.DataFrame) -> pd.DataFrame:
    """
    Table 1: 各パターンで、per-protein の |log fold-error| が小さかったタンパク質数
    行: 2-species allometry / Mean、列: 3パターン
    """
    cols = {}
    for target in ["Monkey", "Rat", "Mouse"]:
        s1, s2, t, label, c_pm, c_pa, c_em, c_ea = PATTERNS[target]
        em = df_out[c_em]
        ea = df_out[c_ea]
        ok = em.notna() & ea.notna()
        cols[label] = {
            "2-species allometry": int((ea[ok] < em[ok]).sum()),
            "Mean": int((em[ok] < ea[ok]).sum()),
        }
        n_tie = int((ea[ok] == em[ok]).sum())
        if n_tie > 0:
            cols[label]["Tie"] = n_tie
    tbl = pd.DataFrame(cols)
    if "Tie" in tbl.index:
        tbl = tbl.fillna(0).astype(int)
    return tbl


def make_table_figure(tbl: pd.DataFrame, dpi: int = 300):
    """Table 1 を論文風（横罫線のみ）の PNG 用 Figure として描画"""
    n_rows = len(tbl.index)
    row_in = 0.42                       # 1行の高さ(inch)
    fig_w = 10.0
    fig_h = 1.5 + row_in * (n_rows + 1)
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)

    y = fig_h - 0.35
    ax.text(0.2, y, "Table 1. Comparison of predictivity between mean and 2-species allometry approaches",
            fontsize=11, fontweight="bold", va="center")

    col_x = [1.4, 4.0, 6.1, 8.3]
    x0, x1 = 0.2, fig_w - 0.2
    y -= 0.55
    ax.plot([x0, x1], [y, y], color="black", lw=1.3)
    y -= row_in / 2
    for x, c in zip(col_x[1:], tbl.columns):
        ax.text(x, y, c, ha="center", va="center", fontsize=10.5)
    y -= row_in / 2
    ax.plot([x0, x1], [y, y], color="black", lw=0.8)
    for r in tbl.index:
        y -= row_in / 2
        ax.text(col_x[0], y, r, ha="center", va="center", fontsize=10.5)
        for x, c in zip(col_x[1:], tbl.columns):
            ax.text(x, y, f"{int(tbl.loc[r, c])}", ha="center", va="center", fontsize=10.5)
        y -= row_in / 2
    ax.plot([x0, x1], [y, y], color="black", lw=1.3)

    y -= 0.3
    ax.text(0.2, y,
            "Each value indicates the number of proteins for which the corresponding approach "
            "(2-species allometry or Mean) yielded lower\nper-protein absolute log$_{10}$ fold-error "
            "in two-species-to-third predictions.",
            fontsize=9, va="top")
    return fig


def fig_to_png_bytes(fig, dpi: int = 300) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return buf.getvalue()


# =========================
# Streamlit UI
# =========================
st.set_page_config(page_title="Scenario prediction: Mean vs Allometry2 (+ Allo3 vs Mean3)", layout="wide")
st.title("半減期予測：2種Mean vs 2点Allometry（3パターン）＋ 3種Mean vs 3点Allometry（追加）")
st.caption("app.py バージョン: 2026-10-27 v3（Figure / Table PNG出力、比較用の基準〔一定値・指数固定 allometry〕対応）")

with st.expander("固定設定（要件/補足）", expanded=True):
    st.markdown(
        f"""
- 入力: Protein / Mouse / Rat / Monkey の半減期（Excel）
- 体重(kg): **Mouse={BW['Mouse']}**, **Rat={BW['Rat']}**, **Monkey={BW['Monkey']}**
- **左側（従来）**: 3パターン（2種→残り1種）
  - Mouse+Rat → Monkey
  - Mouse+Monkey → Rat
  - Rat+Monkey → Mouse
  - 比較: Mean（2種平均） vs Allometry2（2点回帰）
  - 追加: Allometry2の定数 **a,b** もCSVに出力
- **右側（追加 ②）**: 3種Mean（全3種平均） vs 3点Allometry（3点回帰）
  - Mouse/Rat/Monkey をそれぞれ予測（in-sample比較）
  - さらに全体RMSE(log)でも比較結果を出力
- 誤差:
  - 個別: $|\\ln(true)-\\ln(pred)|$
  - 全体: $\\sqrt{{\\mathrm{{mean}}((\\ln(true)-\\ln(pred))^2)}}$
- 注意: 欠損や **0以下** はlogが取れないため、該当の比較は **NA** になります
"""
    )

# 表示ヘルパー（Streamlit バージョン差を吸収）
def show_png(png: bytes):
    """Streamlit のバージョン差（use_container_width / width="stretch"）を吸収して画像表示"""
    try:
        st.image(png, width="stretch")
        return
    except Exception:
        pass
    try:
        st.image(png, use_container_width=True)
        return
    except Exception:
        pass
    st.image(png)


def show_df(df, **kw):
    """st.dataframe のバージョン差を吸収"""
    try:
        st.dataframe(df, width="stretch", **kw)
    except Exception:
        try:
            st.dataframe(df, use_container_width=True, **kw)
        except Exception:
            st.dataframe(df, **kw)


# Sidebar
st.sidebar.header("表示フィルタ")
keep_all3pos = st.sidebar.checkbox(
    "3種すべて正の値（All3Positive）の行だけ表示/出力",
    value=False,
    help="データに0/負/欠損が多い場合、比較可能な行に絞ると集計が見やすくなります。",
)

# Load Excel
st.subheader("1) 入力Excelの読み込み")
uploaded = st.file_uploader("Excelファイル(.xlsx)をアップロード", type=["xlsx"])

df_raw = None
source_label = None

def _read_with_sheet(src, label):
    """シートが複数ある場合は選択できるようにする（既定: CommonList）"""
    xls = pd.ExcelFile(src, engine="openpyxl")
    sheets = xls.sheet_names
    idx = sheets.index(DEFAULT_SHEET) if DEFAULT_SHEET in sheets else 0
    sheet = st.selectbox("解析に使用するシート", sheets, index=idx,
                         help="ProteinList_common.xlsx では、3種すべての値が揃った CommonList を推奨")
    return pd.read_excel(xls, sheet_name=sheet), f"{label}（シート: {sheet}）"


if uploaded is not None:
    df_raw, source_label = _read_with_sheet(uploaded, "アップロードされたファイル")
else:
    p = Path(__file__).resolve().parent / DEFAULT_EXCEL
    if not p.exists():
        p = Path(DEFAULT_EXCEL)
    if p.exists():
        df_raw, source_label = _read_with_sheet(p, f"同一フォルダの既定ファイル: {DEFAULT_EXCEL}")

if df_raw is None:
    st.info("Excelをアップロードするか、app.py と同じフォルダに既定ファイルを置いてください。")
    st.stop()

st.caption(f"読み込み元: {source_label}")
st.write("読み込みデータ（先頭5行）")
show_df(df_raw.head(5))

# Run analysis
st.subheader("2) 解析（3パターン × 2方式比較 + 3種比較）")

try:
    df_out = run_analysis(df_raw)
except Exception as e:
    st.error(f"解析に失敗しました: {e}")
    st.stop()

if keep_all3pos:
    df_out = df_out[df_out["All3Positive"]].copy()

# 比較用の基準（一定値・指数固定 allometry）の列を追加（既定 b = 0.25、CSV にも出力）
df_out = add_baselines(df_out, FIXED_EXPONENT_DEFAULT)

# =========================
# Figure / Table 出力（PNG） ※解析直後に表示
# =========================
st.subheader("3) Figure / Table（PNG出力）")
tab_fig, tab_tbl, tab_base = st.tabs(["Figure（Mean vs Allometry）", "Table 1（優位タンパク質数）", "比較用の基準（補足表）"])

with tab_fig:
    fc1, fc2, fc3, fc4 = st.columns(4)
    with fc1:
        fig_target = st.selectbox(
            "Figure の予測パターン",
            ["Monkey", "Rat", "Mouse"],
            format_func=lambda k: PATTERNS[k][3],
            index=0,
        )
    with fc2:
        ax_lo = st.number_input("散布図 軸下限 (days)", min_value=0.001, value=0.5, step=0.1, format="%.3f")
    with fc3:
        ax_hi = st.number_input("散布図 軸上限 (days)", min_value=0.01, value=200.0, step=10.0)
    with fc4:
        box_mode = st.radio("箱ひげ図の値", ["log10(Pred/Obs)（符号あり）", "|log10(Pred/Obs)|（絶対値）"], index=0)
    fig_dpi = st.slider("出力解像度 (dpi)", min_value=100, max_value=600, value=300, step=50)

    if ax_hi <= ax_lo:
        st.error("軸上限は下限より大きくしてください。")
    else:
        try:
            fig, fig_stats = make_allometry_figure(
                df_out, target=fig_target, lim=(float(ax_lo), float(ax_hi)),
                signed_box=box_mode.startswith("log10"), dpi=fig_dpi,
            )
            fig_png = fig_to_png_bytes(fig, dpi=fig_dpi)
            st.markdown(f"#### Figure: {PATTERNS[fig_target][3]}（Mean vs 2-species allometry）")
            show_png(fig_png)
            st.download_button(
                "Figure をPNGでダウンロード", data=fig_png,
                file_name=f"Allometry_Fig_{fig_target}.png", mime="image/png",
                key="dl_fig",
            )
            st.caption(
                "(a) Mean 予測 vs 観測、(b) Allometry 予測 vs 観測（同一軸範囲）。実線=1:1、破線=2倍、点線=3倍の予測幅。"
                "軸範囲外の予測値は上端/下端に三角で表示。(c) 箱ひげ図（箱=IQR、中央線=中央値、ひげ=1.5×IQR内の最端値、外れ値非表示）。"
                "P値は per-protein の |log10 fold-error| に対する両側 Wilcoxon signed-rank 検定。"
            )
            stats_df = pd.DataFrame([fig_stats]).T.rename(columns={0: "Value"}).astype(str)
            show_df(stats_df)
        except Exception as e:
            st.error(f"Figure の作成に失敗しました: {e}")
            st.exception(e)

with tab_tbl:
    try:
        tbl = make_comparison_table(df_out)
        show_df(tbl)
        tbl_png = fig_to_png_bytes(make_table_figure(tbl, dpi=300), dpi=300)
        show_png(tbl_png)
        tc1, tc2 = st.columns(2)
        with tc1:
            st.download_button(
                "Table をPNGでダウンロード", data=tbl_png,
                file_name="Allometry_Table.png", mime="image/png", key="dl_tbl_png",
            )
        with tc2:
            st.download_button(
                "Table をCSVでダウンロード",
                data=tbl.to_csv(encoding="utf-8-sig").encode("utf-8-sig"),
                file_name="Allometry_Table.csv", mime="text/csv", key="dl_tbl_csv",
            )
    except Exception as e:
        st.error(f"Table の作成に失敗しました: {e}")
        st.exception(e)


with tab_base:
    st.markdown(
        "予測法の成績を、ターゲット種のデータを使わない2つの基準と同じタンパク質集合で比較します。  \n"
        "**① 一定値の予測**: 全タンパク質に同じ値（2種平均予測の、全タンパク質での中央値）を当てる。"
        "タンパク質ごとの情報を使わずにどれだけ当たるかの物差し。  \n"
        "**② 指数固定の allometry**: source 2種の半減期をそれぞれ (BW_target / BW_source)^b で補正し、その幾何平均を予測値とする。"
        "b = 0.25 は生理学的時間のスケーリング（t ∝ BW^0.25; Boxenbaum 1982）。"
    )
    bc1, bc2 = st.columns(2)
    with bc1:
        base_target = st.selectbox("予測パターン", ["Monkey", "Rat", "Mouse"],
                                   format_func=lambda k: PATTERNS[k][3], index=0, key="base_target")
    with bc2:
        base_exp = st.number_input("指数 b（固定）", value=float(FIXED_EXPONENT_DEFAULT), step=0.05, format="%.2f",
                                   help="既定 0.25（生理学的時間のスケーリング）")
    try:
        df_b = df_out if abs(base_exp - FIXED_EXPONENT_DEFAULT) < 1e-12 else add_baselines(df_out, base_exp)
        btbl = make_baseline_table(df_b, base_target)
        const_v = df_b.attrs.get(f"Const_{base_target}", np.nan)
        st.caption(f"一定値 = {const_v:.3f} days（{PATTERNS[base_target][3]} の2種平均予測の中央値）、指数 b = {base_exp:g}")
        show_df(btbl.round(3))
        bpng = fig_to_png_bytes(make_baseline_table_figure(btbl, const_v, base_exp, dpi=300), dpi=300)
        show_png(bpng)
        bd1, bd2 = st.columns(2)
        with bd1:
            st.download_button("補足表をPNGでダウンロード", data=bpng,
                               file_name=f"Baseline_Table_{base_target}.png", mime="image/png", key="dl_base_png")
        with bd2:
            st.download_button("補足表をCSVでダウンロード",
                               data=btbl.to_csv(encoding="utf-8-sig").encode("utf-8-sig"),
                               file_name=f"Baseline_Table_{base_target}.csv", mime="text/csv", key="dl_base_csv")
    except Exception as e:
        st.error(f"比較用の基準の計算に失敗しました: {e}")
        st.exception(e)


# Summary: Best method counts for each 2-species pattern
st.subheader("参考: 詳細集計")
st.markdown("### 2種→1種（3パターン）: Mean vs Allometry2 のベスト件数")
c1, c2, c3 = st.columns(3)

with c1:
    st.write("Mouse+Rat → Monkey")
    show_df(
        df_out["BestMonkey_MeanVsAllo2"].value_counts(dropna=False).rename("Count").to_frame(),
    )

with c2:
    st.write("Mouse+Monkey → Rat")
    show_df(
        df_out["BestRat_MeanVsAllo2"].value_counts(dropna=False).rename("Count").to_frame(),
    )

with c3:
    st.write("Rat+Monkey → Mouse")
    show_df(
        df_out["BestMouse_MeanVsAllo2"].value_counts(dropna=False).rename("Count").to_frame(),
    )

# Summary: 3-species comparison
st.markdown("### 3種（追加②）: Mean3 vs Allometry3 の比較（件数）")
c4, c5, c6, c7 = st.columns(4)

with c4:
    st.write("Mouse（3種モデル）")
    show_df(
        df_out["BestMouse_Mean3VsAllo3"].value_counts(dropna=False).rename("Count").to_frame(),
    )

with c5:
    st.write("Rat（3種モデル）")
    show_df(
        df_out["BestRat_Mean3VsAllo3"].value_counts(dropna=False).rename("Count").to_frame(),
    )

with c6:
    st.write("Monkey（3種モデル）")
    show_df(
        df_out["BestMonkey_Mean3VsAllo3"].value_counts(dropna=False).rename("Count").to_frame(),
    )

with c7:
    st.write("Overall（RMSE log）")
    show_df(
        df_out["BestOverall_Mean3VsAllo3"].value_counts(dropna=False).rename("Count").to_frame(),
    )

# Table
st.markdown("### 行ごとの結果（全列）")
show_df(df_out, height=650)

# Export CSV
st.subheader("4) CSVエクスポート")
csv_bytes = df_out.to_csv(index=False, encoding="utf-8-sig").encode("utf-8-sig")

st.download_button(
    label="結果CSVをダウンロード",
    data=csv_bytes,
    file_name="results.csv",
    mime="text/csv",
)

st.caption("※ Excelで文字化けしにくいよう UTF-8-SIG で出力します。")

