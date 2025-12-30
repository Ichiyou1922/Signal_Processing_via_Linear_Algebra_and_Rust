import sys
import pandas as pd
import matplotlib.pyplot as plt

# 標準入力からデータを読み込み
try:
    df = pd.read_csv(sys.stdin)
except Exception as e:
    print("Error reading CSV:", e)
    sys.exit(1)

# プロット設定
plt.figure(figsize=(10, 6))

# 矩形窓 (Rectangular) - 青色・破線
plt.plot(df['k'], df['mag_rect'], label='Rectangular (No Window)', 
         marker='o', linestyle='--', color='blue', alpha=0.6)

# ハニング窓 (Hanning) - 赤色・実線
plt.plot(df['k'], df['mag_hann'], label='Hanning Window', 
         marker='s', linestyle='-', color='red', linewidth=2)

# グラフ装飾
plt.title('Spectrum Leakage Analysis: Rectangular vs Hanning')
plt.xlabel('Frequency Index k')
plt.ylabel('Magnitude (Linear)')
plt.grid(True)
plt.legend()
plt.xticks(df['k'][::2]) # 目盛りを適度に間引く

# 保存または表示
plt.savefig('spectrum_analysis.png')
print("Graph saved to 'spectrum_analysis.png'")
# plt.show()
