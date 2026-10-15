// Keep low-level exception prefixes out of player-facing feedback.
export function playerError(message) {
  const text=String(message||'').replace(/^(?:std::(?:runtime_error|invalid_argument):\s*)+/,'');
  if (/Failed to fetch|NetworkError|Load failed|network.*(?:failed|error)/i.test(text))
    return '资源下载失败，请检查网络连接，然后再次点击“保存并进入赛道”重试。';
  if (/QuotaExceededError|Storage quota exceeded/i.test(text))
    return '浏览器暂存空间不足，请点击“保存工程”下载备份。';
  return text;
}
